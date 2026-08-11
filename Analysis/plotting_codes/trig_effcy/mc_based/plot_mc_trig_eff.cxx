// plot_mc_trig_eff.cxx
// Step-1/2/3 MC trigger-efficiency plots (mc_trigger_efficiency.md §3.1-§3.3).
//
//   Step 1 (§3.1): MC direct P[mu4-matched | offline reco muon] vs data tag-and-probe
//                  P[2mu4 | mu4-tag, DeltaR>0.8] (the data eps^nc estimate), per charge.
//   Step 2 (§3.2): DeltaR-binned MC singles efficiency (factorization cross-check),
//                  3 DeltaR bins overlaid + Step-1 inclusive reference.
//   Step 3 (§3.3): eps_dR(dR) = inverse-weighted num / unweighted denom; plateau = weighted
//                  mean over dR in [1,4].
//   Step 4 (§3.4): eps_single(dR), the leg-level analog.
//
// ERROR BARS on the Step-3/Step-4 ratios (round 7): the numerator is a re-weighted SUBSET of
// the denominator, so TH1::Divide's independent propagation over-states them by
// sqrt((1+eps)/(1-eps)) = 1.5-3.2x. They now use the conditional (binomial-correct) form of
// SetConditionalRatioErrors() below, fed by the errA/errB (+ Step-4 covP/covQ) histograms
// booked in FillMCTrigEffHists.cxx. Central values are unchanged.
//
// One macro for both samples:
//   plot_mc_trig_eff("pp")      : Pythia8 pp24 fullsim  vs pp24 data       -> eps_dR^2mu4
//   plot_mc_trig_eff("overlay") : HIJING overlay PbPb23 vs PbPb23 data 0-5% (D2) -> eps_dR^cross
//
// WP config (registry: Analysis/docs/muon_wp_registry.md): use_tight_wp default TRUE (Tight
// nominal, unsuffixed inputs); false selects the _medium_wp MC inputs AND the WP-matched
// _medium_wp DATA tag-and-probe file (+ label text). MC and data are never compared across
// working points.
//
// Data sign convention (evidence in docs/tracking/_sub_mctrig_plots.md F3):
//   sign1 = mu+ (RDFBasedHistFillingPP.cxx:168 Filter("charge2nd > 0");
//   SingleMuEffcyPtTurnOnFitter.cxx:153-154), sign2 = mu-.
//
// Compile/run (ACLiC, from this directory):
//   root -l -b -q 'plot_mc_trig_eff.cxx+("pp")'
//   root -l -b -q 'plot_mc_trig_eff.cxx+("overlay")'

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TArrow.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TNamed.h>
#include <TROOT.h>

#include <TH3D.h>
#include <TGraphErrors.h>
#include <TPad.h>

// Sample identity (input dir, file label, plot root, headline, eps symbol, FULL-vs-TEST flag)
// and the names of the plateau / fit ROOT files. Shared with fit_dr_corrections.cxx and
// plot_dr_correction_fits.cxx so the measurement and the fit can never drift apart.
#include "dr_correction_sample_cfg.h"
#include "../../../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../../../Utilities/proj_range_to_suffix.cxx"
// The round-7 conditional ratio error + the (pair pT, pair eta) cell projection, shared with
// the fit stage so both use the identical error definition.
#include "dr_correction_ratio.h"
// The Step-1 sanity-check pT-match threshold, shared with FillMCTrigEffHists.cxx (which APPLIES
// it) so the value DRAWN on the canvas is the value that was cut on.
#include "../../../Utilities/MCTrigEffSanityCfg.h"
#include "../../../Utilities/MCTrigEffPlateauWindow.h"
#include "../../../Utilities/MCTrigEffPairPtBinning.h"

#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// ---------------------------------------------------------------- helpers

template <typename T>
T* GetObj(TFile* f, const std::string& name)
{
    T* obj = dynamic_cast<T*>(f->Get(name.c_str()));
    if (!obj)
        throw std::runtime_error("plot_mc_trig_eff: missing object '" + name +
                                 "' in file " + f->GetName());
    return obj;
}

// Soft variant: returns nullptr instead of throwing. Used ONLY for the round-9 per-sign
// histograms, which a sample filled before round 9 simply does not have -- a missing per-sign
// series must degrade to "no sign-separated results for this sample", never kill the whole
// plot set. Anything the nominal results depend on still goes through GetObj.
template <typename T>
T* GetObjOrNull(TFile* f, const std::string& name)
{
    return dynamic_cast<T*>(f->Get(name.c_str()));
}

TFile* OpenFile(const std::string& path)
{
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie())
        throw std::runtime_error("plot_mc_trig_eff: cannot open " + path);
    return f;
}

// Bayesian efficiency graph (efficiencies bounded in [0,1]; num subset of denom).
TGraphAsymmErrors* BayesEff(TH1* num, TH1* denom)
{
    auto* g = new TGraphAsymmErrors();
    g->Divide(num, denom, "cl=0.683 b(1,1) mode");
    return g;
}

void StyleGraph(TGraphAsymmErrors* g, Color_t col, Style_t marker,
                double msize = 0.9, Width_t lw = 2, Style_t ls = 1)
{
    g->SetMarkerColor(col);
    g->SetLineColor(col);
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(msize);
    g->SetLineWidth(lw);
    g->SetLineStyle(ls);
}

// suppress_xlabels: for the MAIN pad of a split (main+ratio) canvas -- the ratio pad below
// carries the x axis, and leaving the main pad's x labels on renders them chopped in half
// against its 2% bottom margin.
TH1* DrawEffFrame(double xlo, double xhi, const std::string& xtitle,
                  double ylo = 0.0, double yhi = 1.1,
                  const std::string& ytitle = "efficiency",
                  bool suppress_xlabels = false)
{
    TH1* frame = gPad->DrawFrame(xlo, ylo, xhi, yhi);
    if (gPad->GetLogx()) { frame->GetXaxis()->SetMoreLogLabels(); frame->GetXaxis()->SetNoExponent(); }
    frame->GetXaxis()->SetTitle(xtitle.c_str());
    frame->GetYaxis()->SetTitle(ytitle.c_str());
    frame->GetXaxis()->SetTitleSize(suppress_xlabels ? 0.0 : 0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetLabelSize(suppress_xlabels ? 0.0 : 0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    return frame;
}

// ---- ratio panels (review criterion R3) --------------------------------------------
// Step 1 overlays MC on a DATA reference and Step 2 overlays three DeltaR series that are
// EXPECTED TO AGREE (the §3.2 factorization cross-check) -- in both cases the ratio IS the
// measurement, so every such canvas carries a bottom ratio pad.

// Split the CURRENT pad into a main pad (top) and a ratio pad (bottom). Returns {main, ratio}.
// Both pads are drawn; caller cd()s into them.
std::pair<TPad*, TPad*> SplitPadForRatio(const std::string& tag, bool logx,
                                         double split = 0.32)
{
    TPad* host = (TPad*)gPad;
    host->cd();
    auto* pmain = new TPad(("pm_" + tag).c_str(), "", 0., split, 1., 1.);
    auto* prat  = new TPad(("pr_" + tag).c_str(), "", 0., 0.,    1., split);
    pmain->SetTopMargin(0.06);
    pmain->SetBottomMargin(0.02);
    pmain->SetLeftMargin(0.14);
    prat->SetTopMargin(0.03);
    prat->SetBottomMargin(0.35);
    prat->SetLeftMargin(0.14);
    if (logx) { pmain->SetLogx(); prat->SetLogx(); }
    pmain->Draw();
    prat->Draw();
    return {pmain, prat};
}

// Frame for a ratio pad: labels/titles scaled up by 1/split so they match the main pad.
TH1* DrawRatioFrame(double xlo, double xhi, const std::string& xtitle,
                    const std::string& ytitle, double ylo, double yhi,
                    double split = 0.32)
{
    const double s = 1.0 / split;
    TH1* fr = gPad->DrawFrame(xlo, ylo, xhi, yhi);
    if (gPad->GetLogx()) { fr->GetXaxis()->SetMoreLogLabels(); fr->GetXaxis()->SetNoExponent(); }
    fr->GetXaxis()->SetTitle(xtitle.c_str());
    fr->GetYaxis()->SetTitle(ytitle.c_str());
    fr->GetXaxis()->SetTitleSize(0.05 * s);
    fr->GetYaxis()->SetTitleSize(0.05 * s);
    fr->GetXaxis()->SetLabelSize(0.045 * s);
    fr->GetYaxis()->SetLabelSize(0.045 * s);
    fr->GetYaxis()->SetTitleOffset(0.42);
    fr->GetXaxis()->SetTitleOffset(1.05);
    fr->GetYaxis()->SetNdivisions(505);
    auto* l = new TLine(xlo, 1.0, xhi, 1.0);
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->Draw("same");
    return fr;
}

// The conditional (binomial-correct) ratio error of round 7 now lives in
// dr_correction_ratio.h (included above), so the fit stage judges chi2/ndf with exactly the
// same errors these plots show.

// Point-by-point ratio of two efficiency graphs.
//
// MATCH BY X, NOT BY INDEX. TGraphAsymmErrors::Divide SKIPS bins with an empty denominator,
// so the numerator and denominator graphs do NOT share an index->bin mapping: one skipped
// bin shifts every later index. An index-paired loop then either divides mismatched bins or
// (with an x-guard) silently DROPS every point after the first skip -- which cost 12 real
// ratio points per working point in pp. Pairing on the x value is immune to this.
//
// The ONLY point excluded is one whose DENOMINATOR is zero: the ratio is then genuinely
// undefined. (The degenerate zero-width 8.0 GeV pT bin never reaches this guard --
// TGraphAsymmErrors::Divide already drops it upstream.) A zero NUMERATOR with a nonzero
// denominator is a real measurement of ratio 0
// and IS emitted, so the caller's MarkOffScale() gives it a down-arrow instead of letting it
// vanish. (It does fire: e.g. pp mu-, 0.2<=dR<1.0 at q.eta=2.39.)
//
// Errors: relative errors in quadrature, symmetrized -- adequate for a ratio pad.
TGraphAsymmErrors* DivideGraphClean(TGraphAsymmErrors* gn, TGraphAsymmErrors* gd)
{
    auto* gr = new TGraphAsymmErrors();
    int k = 0;
    for (int i = 0; i < gn->GetN(); ++i) {
        double xn, yn;
        gn->GetPoint(i, xn, yn);

        // find the denominator point at the same x
        int j = -1;
        for (int m = 0; m < gd->GetN(); ++m) {
            double xd, yd;
            gd->GetPoint(m, xd, yd);
            if (std::fabs(xn - xd) <= 1e-6 * std::max(1.0, std::fabs(xd))) { j = m; break; }
        }
        if (j < 0) continue;                       // no denominator bin at this x

        double xd, yd;
        gd->GetPoint(j, xd, yd);
        if (yd <= 0.) continue;                    // ratio undefined

        const double en = 0.5 * (gn->GetErrorYhigh(i) + gn->GetErrorYlow(i));
        const double ed = 0.5 * (gd->GetErrorYhigh(j) + gd->GetErrorYlow(j));
        const double r  = yn / yd;
        // relative error of the numerator is undefined at yn == 0; use the denominator's
        // alone (the point is at r = 0 and only needs to be visible, not precise).
        const double rel_n = (yn > 0.) ? (en / yn) : 0.;
        const double er = (yn > 0.) ? r * std::sqrt(rel_n * rel_n + (ed / yd) * (ed / yd))
                                    : en / yd;
        gr->SetPoint(k, xn, r);
        gr->SetPointError(k, gn->GetErrorXlow(i), gn->GetErrorXhigh(i), er, er);
        ++k;
    }
    return gr;
}

// A zoomed ratio frame hides points that fall outside it. Mark every such point with an
// arrow at the frame edge, so a clipped point can never be mistaken for a missing one --
// the same rule the Step-3 slice panel follows. (Points whose CENTRAL value is off-scale;
// an error bar running past the frame with the central value inside is fine.)
void MarkOffScale(TGraphAsymmErrors* g, double ylo, double yhi, Color_t col)
{
    for (int i = 0; i < g->GetN(); ++i) {
        double x, y;
        g->GetPoint(i, x, y);
        if (y <= yhi && y >= ylo) continue;
        const bool up  = (y > yhi);
        const double y0 = up ? ylo + 0.88 * (yhi - ylo) : ylo + 0.12 * (yhi - ylo);
        const double y1 = up ? ylo + 0.98 * (yhi - ylo) : ylo + 0.02 * (yhi - ylo);
        auto* ar = new TArrow(x, y0, x, y1, 0.010, "|>");
        ar->SetLineColor(col);
        ar->SetFillColor(col);
        ar->SetLineWidth(2);
        ar->Draw();
    }
}

void DrawUnityLine(double xlo, double xhi)
{
    auto* l = new TLine(xlo, 1.0, xhi, 1.0);
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->Draw("same");
}

void SaveCanvas(TCanvas& c, const std::string& path)
{
    c.SaveAs(path.c_str());
    std::cout << "  wrote " << path << std::endl;
}

// ---------------------------------------------------------------- config

struct SampleCfg {
    std::string mc_dir;
    std::string mc_label;      // file label
    std::string data_file;
    std::string ctr;           // data ctr token ("" or "_ctr0_5")
    std::string out_base;      // plots root for this sample
    std::string sample_text;   // canvas headline (without WP)
    std::string data_text;     // data legend entry
    std::string eps_dr_text;   // eps_DR symbol for step 3
    bool step2_coarse;         // rebin the Step-2 DeltaR-comparison panels (see below)

    // OPTIONAL second MC sample, overlaid on the Step-1 q.eta-binned panels ONLY, as
    // points with NO fit curve. Used by "noovl" to put the HIJING overlay (r17618) beside
    // r17663 -- that MC-vs-MC comparison IS the deliverable of the r17663 study (R10).
    // Empty cmp_fit_file => no third series (pp and overlay keep their two-series panels
    // exactly as reviewed).
    std::string cmp_fit_file;  // file holding g_mc_pt_vs_q_eta_<charge>_<qeta> for the comparison MC
    std::string cmp_text;      // legend entry

    // OPTIONAL: on the Step-1 q.eta-binned panels ONLY, replace the black DATA tag-and-probe
    // series with a black MC series (points, no fit) from this fit file's
    // g_mc_pt_vs_q_eta_<charge>_<qeta> graphs. Used by "noovl" to make the panel a pure
    // MC-vs-MC-vs-MC comparison (r17663 vs pp24-conditions vs HIJING-overlay). When set, the
    // q.eta-panel ratio pad becomes <this sample MC> / <black MC> and the headline/legend say
    // so. Empty => the black series is the data T&P reference, exactly as pp/overlay use it.
    std::string qeta_black_mc_file;  // pp24 fullsim MC fit file for the black series
    std::string qeta_black_mc_text;  // legend entry for the black MC series
};

// The DATA reference must be at the SAME working point as the MC (§3.0(d)): the data
// tag-and-probe hist file carries the WP in its name (`_medium_wp` when the data RDF ran
// with isTight=false; unsuffixed = Tight nominal). Comparing Medium MC against the Tight
// data file -- which is what a hardcoded path would silently do -- is meaningless.
SampleCfg MakeCfg(const std::string& sample, bool use_tight_wp)
{
    const std::string data_wp = use_tight_wp ? "" : "_medium_wp";

    // Sample IDENTITY (mc_dir / mc_label / out_base / sample_text / eps_dr_text) comes from the
    // shared table in dr_correction_sample_cfg.h -- the same table the fit stage reads, so a
    // path or label can never drift between the step that MEASURES the plateau and the step
    // that CONSUMES it. Everything below is plot-macro-specific and stays here.
    const DrCorrSample id = GetDrCorrSample(sample, use_tight_wp);

    SampleCfg c;
    c.mc_dir      = id.mc_dir;
    c.mc_label    = id.mc_label;
    c.out_base    = id.out_base;
    c.sample_text = id.sample_text;
    c.eps_dr_text = id.eps_dr_text;

    if (sample == "pp") {
        c.data_file   = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
                        "histograms_real_pairs_pp_2024_single_mu4_coarse_q_eta_bin_qeta_fid"
                        + data_wp + ".root";
        c.ctr         = "";
        c.data_text   = "pp 2024 data";
        c.step2_coarse = false;  // pp has ~4x the pair statistics: the fine axes are readable
    } else if (sample == "pp_full") {
        // pp24 FULL sample. Physics identical to "pp" -- SAME pp24 data reference, SAME
        // 2mu4 product weighting, SAME plot LOCATION (this is the canonical pp trig-eff
        // deliverable, which the full sample now supersedes -- back up the TEST-sample plots
        // first, done by the pipeline / the run wrapper). ONLY the MC inputs differ:
        // mc_dir = full-sample dir, mc_label = "pp24_full" (reads the _full intermediate hists,
        // so the hists/fits never clobber the TEST ones).
        c.data_file   = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
                        "histograms_real_pairs_pp_2024_single_mu4_coarse_q_eta_bin_qeta_fid"
                        + data_wp + ".root";
        c.ctr         = "";
        c.data_text   = "pp 2024 data";
        c.step2_coarse = false;  // the full sample has far MORE pair statistics than the test
    } else if (sample == "overlay") {
        c.data_file   = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/"
                        "histograms_real_pairs_pbpb_2023_single_mu4_coarse_q_eta_bin_qeta_fid"
                        + data_wp + ".root";
        c.ctr         = "_ctr0_5";   // D2: overlay compares ONLY to PbPb23 data 0-5%
        c.data_text   = "Pb+Pb 2023 data, 0-5%";
        // The overlay pair sample (0-5% only) is ~4x thinner than pp; on the native
        // 41-bin pT / 184-bin q.eta axes the three DeltaR series are an unreadable
        // error-bar forest and the comparison the panel exists for cannot be made.
        c.step2_coarse = true;
    } else if (sample == "noovl") {
        // r17663 NO-OVERLAY diagnostic (R8, round 4). Its own plot root -- the pp24 and
        // overlay plot sets are NOT touched.
        //
        // DATA REFERENCE = pp24 data, deliberately: r17663 simulates pp COLLISIONS with no
        // overlaid event, so the pp data turn-on is the like-for-like measurement (the same
        // reference the pp24 fullsim is validated against, which is what makes the two MC
        // curves directly comparable). Its RECO CONDITIONS are PbPb23-like, which is exactly
        // the variable under test -- so the data curve here is CONTEXT, not the deliverable:
        // the deliverable is the MC-vs-MC comparison of the forward bin (R8 outcome tree).
        c.data_file   = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
                        "histograms_real_pairs_pp_2024_single_mu4_coarse_q_eta_bin_qeta_fid"
                        + data_wp + ".root";
        c.ctr         = "";          // no centrality: there is no overlaid event
        c.data_text   = "pp 2024 data";
        // Single 10k-event slice: statistics are thinner than pp24 -> coarse Step-2 axes,
        // for the same readability reason as the overlay.
        c.step2_coarse = true;
        // Third series on the Step-1 q.eta panels: the HIJING overlay (r17618). r17663 vs
        // r17618 is the MC-vs-MC comparison the sample exists for -- r17663 differs from
        // r17618 ONLY by the absence of HIJING, so the two curves together separate
        // "occupancy" from "r16578 configuration" (R10). Points only, NO fit: the overlay's
        // own turn-on fit belongs to its own plot set (and its fit mode is fermi+log, not
        // the erf+log used here) -- drawing it would invite a fit-quality reading of a
        // curve that is here purely as a reference sample.
        c.cmp_fit_file = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/"
                         "single_mu_effcy_pT_fit_mc" + std::string(use_tight_wp ? "" : "_medium_wp")
                       + ".root";
        c.cmp_text     = "MC, Pb+Pb 2023 conditions, with HIJING overlay";
        // Black series on the q.eta panels = pp24-CONDITIONS fullsim MC (FULL sample), NOT
        // pp data: the r17663 study is an MC-vs-MC comparison of reco-tag CONFIGURATIONS, so
        // all three curves are MC (red = Pb+Pb23 cond. no overlay; black = pp24 cond.;
        // blue = Pb+Pb23 cond. with HIJING overlay). Full sample chosen for statistics.
        // FitMCSinglesEffcy writes the unqualified basename, distinguished here only by the
        // full_sample directory (Remaining Work 3c).
        c.qeta_black_mc_file = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/"
                               "single_mu_effcy_pT_fit_mc"
                             + std::string(use_tight_wp ? "" : "_medium_wp") + ".root";
        c.qeta_black_mc_text = "MC, pp 2024 conditions";
    } else {
        throw std::runtime_error("plot_mc_trig_eff: sample must be 'pp', 'pp_full', 'overlay' or "
                                 "'noovl', got " + sample);
    }
    return c;
}

// charge <-> data-sign mapping (F3): sign1 = mu+, sign2 = mu-
const std::vector<std::string> kCharges     = {"muplus", "muminus"};
const std::vector<std::string> kDataSigns   = {"sign1", "sign2"};
const std::vector<std::string> kChargeTex   = {"#mu^{+}", "#mu^{-}"};

// q.eta bins -- DERIVED from CommonEffcyConfig, never retyped. Two hand-maintained copies used
// to live here (a suffix list and a numeric list) and they are exactly the kind of duplicate the
// binning rule exists to prevent: they had to be edited in lockstep with the fitter's own third
// copy every time the binning moved. Since round 8 the NOMINAL binning is the CONTIGUOUS COARSE
// one (gaps included), so there are no unfitted holes.
const std::vector<std::pair<double,double>>& kQEtaRangeF = [] () -> const std::vector<std::pair<double,double>>& {
    static std::vector<std::pair<double,double>> v;
    static const CommonEffcyConfig cfg{};
    if (v.empty()) for (const auto& r : cfg.q_eta_proj_ranges_coarse_incl_gap) v.emplace_back(r.first, r.second);
    return v;
}();
const std::vector<std::pair<double,double>>& kQEtaRange = kQEtaRangeF;
const std::vector<std::string> kQEtaSuffix = [] {
    static const CommonEffcyConfig cfg{};
    std::vector<std::string> v;
    for (const auto& r : cfg.q_eta_proj_ranges_coarse_incl_gap) v.push_back(pairToSuffix(r));
    return v;
}();

// ---- Step-2 coarse binning (SampleCfg::step2_coarse) -------------------------------
// Step 2 asks ONE question: at FIXED (pT, q.eta), do the three DeltaR series agree?
// That needs bins with enough entries to separate the series, not the fine resolution of
// the Step-1 turn-on fit. So merge the native bins into coarse ones that still sample the
// turn-on and collapse the sparse high-pT tail.
//
// TH1::Rebin REQUIRES every new edge to coincide with an existing one. A rounded literal
// (5.1 for the native 5.09824) does NOT coincide: ROOT then prints "Bin edge ... does not
// match ... Result can be inconsistent" and groups by bin CENTRE, which can split a native
// bin. So the requested values below are only TARGETS -- SnapToAxis() replaces each with
// the nearest actual edge of the histogram being rebinned, and THROWS if none is close.
// Nothing downstream ever sees a rounded edge.
const std::vector<double> kStep2CoarsePtTarget =
    {4.0, 5.1, 6.5, 8.0, 10.8, 16.2, 26.8, 60.0};
std::vector<double> Step2CoarseQEtaTarget()   // uniform 0.2 (exact edges of the 184-bin axis)
{
    std::vector<double> e;
    for (int i = 0; i <= 24; ++i) e.push_back(-2.4 + 0.2 * i);
    return e;
}

// Snap each target to the nearest EXACT edge of h's x-axis. Throws if the nearest edge is
// further than tol_rel (relative), i.e. if the requested binning is not achievable -- that
// is a coding error, never something to paper over with an inconsistent rebin.
std::vector<double> SnapToAxis(const TH1* h, const std::vector<double>& targets,
                               double tol_rel = 0.05)
{
    const TAxis* ax = h->GetXaxis();
    std::vector<double> out;
    for (double t : targets) {
        double best = 0.;
        double bestd = 1e300;
        for (int i = 1; i <= ax->GetNbins() + 1; ++i) {
            const double e = ax->GetBinLowEdge(i);
            const double d = std::fabs(e - t);
            if (d < bestd) { bestd = d; best = e; }
        }
        const double scale = std::max(1.0, std::fabs(t));
        if (bestd > tol_rel * scale)
            throw std::runtime_error(Form("SnapToAxis: no axis edge near %.4g on '%s' "
                                          "(nearest %.4g)", t, h->GetName(), best));
        if (!out.empty() && best <= out.back())
            throw std::runtime_error(Form("SnapToAxis: edges collapsed at %.4g on '%s' "
                                          "(coarse binning finer than the native axis)",
                                          t, h->GetName()));
        out.push_back(best);
    }
    return out;
}

// Rebin to the given TARGET edges (snapped to the native axis first). Returns a clone; the
// input file histogram is never modified. Rebinning num and denom identically preserves
// num <= denom, so the Bayes divide stays valid.
TH1* RebinTo(TH1* h, const std::vector<double>& targets, const std::string& name)
{
    const std::vector<double> edges = SnapToAxis(h, targets);
    return h->Rebin(static_cast<int>(edges.size()) - 1, name.c_str(), edges.data());
}

// DeltaR bins of Step 2
const std::vector<std::string> kDrSuffix = {"dr0_0_2", "dr0_2_1_0", "dr1_0_inf"};
const std::vector<std::string> kDrTex    = {"#DeltaR < 0.2", "0.2 #leq #DeltaR < 1.0",
                                            "#DeltaR #geq 1.0"};
const std::vector<Color_t>     kDrColor  = {kRed + 1, kBlue + 1, kGreen + 2};
const std::vector<Style_t>     kDrMarker = {20, 21, 22};

const Color_t kMCColor   = kRed + 1;
const Color_t kDataColor = kBlack;
// optional comparison MC (SampleCfg::cmp_fit_file) -- blue, distinct from both
const Color_t kCmpColor  = kBlue + 1;

// Plateau window + its systematic variation: SINGLE SOURCE OF TRUTH in
// Analysis/Utilities/MCTrigEffPlateauWindow.h (that header carries the full physics
// justification for both edges). Never retype the edges here -- three separate copies of them
// silently kept the old values when the window changed on 2026-08-04.
const double kPlateauLo     = MCTrigEffPlateau::kLo;
const double kPlateauHi     = MCTrigEffPlateau::kHi;
const double kPlateauSystLo = MCTrigEffPlateau::kSystLo;
const double kPlateauSystHi = MCTrigEffPlateau::kSystHi;

// The overlay headline ("HIJING overlay Pb+Pb23 cond., 0-5%, Tight muons, #mu^{+}") is ~2x
// the length of the pp one and clipped at the pad edge -- taking the charge with it, which
// is the one thing the reader cannot afford to lose. Scale the text down for long strings
// instead of hard-coding a size per call site.
// Count RENDERED glyphs, not raw characters: "#mu^{+}" is 7 characters but one glyph, so
// scaling on text.size() over-shrinks any headline carrying LaTeX markup.
size_t GlyphLength(const std::string& s)
{
    size_t n = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '{' || s[i] == '}' || s[i] == '^' || s[i] == '_') continue;
        if (s[i] == '#') {                       // a control word: #mu, #DeltaR, ... = 1 glyph
            ++i;
            while (i < s.size() && std::isalpha(static_cast<unsigned char>(s[i]))) ++i;
            --i;
        }
        ++n;
    }
    return n;
}

void DrawHeadline(const std::string& text, double x = 0.12, double y = 0.955,
                  double size = 0.038)
{
    constexpr size_t kFitsAt = 34;      // glyphs that fit at the nominal size
    const size_t len = GlyphLength(text);
    if (len > kFitsAt)
        size *= static_cast<double>(kFitsAt) / static_cast<double>(len);
    size = std::max(size, 0.030);       // floor: below this the charge superscript is illegible
    TLatex tl;
    tl.SetNDC();
    tl.SetTextSize(size);
    tl.SetTextFont(42);
    tl.DrawLatex(x, y, text.c_str());
}

// weighted mean of ratio-hist bins with centers in [xlo,xhi]; returns {mean, err}
std::pair<double,double> PlateauWeightedMean(TH1* ratio, double xlo, double xhi)
{
    double sumw = 0., sumwv = 0.;
    for (int i = 1; i <= ratio->GetNbinsX(); ++i) {
        const double xc = ratio->GetBinCenter(i);
        if (xc < xlo || xc > xhi) continue;
        const double v = ratio->GetBinContent(i);
        const double e = ratio->GetBinError(i);
        if (e <= 0. || v == 0.) continue;
        const double w = 1. / (e * e);
        sumw  += w;
        sumwv += w * v;
    }
    if (sumw <= 0.) throw std::runtime_error("PlateauWeightedMean: no usable bins in window");
    return {sumwv / sumw, std::sqrt(1. / sumw)};
}

// ================================================================================
// MACHINE-READABLE PLATEAU OUTPUT (round 7)
//
// The large-dR plateau of eps_dR (Step 3) and eps_single (Step 4) is MEASURED here, per
// (pair pT, pair eta) cell, and CONSUMED by the fit stage (fit_dr_corrections.cxx), which
// divides each cell's dR curve by it before fitting. It therefore has to travel as DATA:
// written to a ROOT file by the step that measures it, read back by the step that uses it.
// The human-readable .txt tables beside the plots stay, but no code may ever parse them --
// a number retyped out of a table (or out of a tracking .md) is a silent-wrong-result waiting
// to happen the next time the chain is re-run.
//
// The 2D maps clone their axes from the source TH3D, so the consumer's (pair pT, pair eta)
// cell indexing is identical to the producer's by construction, not by convention.
// ================================================================================
// `syst` = |mean - mean(retired [1,4] window)|, the plateau-window normalization systematic.
// It is < 0 when the alternative window has no usable dR bin, i.e. "not evaluable" -- the
// consumer must not silently read that as zero uncertainty.
struct PlateauCell { double mean, err, rms; int nb; double syst = -1.; };

// The per-cell large-dR plateau of a dR-correction ratio histogram, measured in the nominal
// window and in the retired [1,4] one (same histogram, so their difference is purely the window).
// Extracted from the Step-3 table block (round 9) so the sign-integrated, same-sign and
// opposite-sign series are all measured by ONE piece of code -- three copies of a weighted mean
// is exactly how two plateau definitions would start to drift apart.
PlateauCell PlateauFromRatio(const TH1D* r)
{
    auto window = [&](double lo, double hi) -> PlateauCell {
        double sw = 0, swv = 0;
        std::vector<std::pair<double,double>> vw;  // (value, weight)
        for (int i = 1; i <= r->GetNbinsX(); ++i) {
            const double xc = r->GetBinCenter(i);
            if (xc < lo || xc > hi) continue;
            const double v = r->GetBinContent(i), e = r->GetBinError(i);
            if (e <= 0. || v == 0.) continue;
            const double w = 1. / (e * e); sw += w; swv += w * v; vw.push_back({v, w});
        }
        if (sw <= 0.) return PlateauCell{ -1, -1, -1, 0 };
        const double mean = swv / sw, err = std::sqrt(1. / sw);
        double swd = 0;
        for (auto& q : vw) swd += q.second * (q.first - mean) * (q.first - mean);
        const double rms = std::sqrt(swd / sw);  // weighted RMS scatter about the mean
        return PlateauCell{ mean, err, rms, (int)vw.size() };
    };
    PlateauCell nom = window(kPlateauLo, kPlateauHi);
    const PlateauCell alt = window(kPlateauSystLo, kPlateauSystHi);
    if (nom.nb > 0 && alt.nb > 0) nom.syst = std::fabs(nom.mean - alt.mean);
    return nom;
}

TH2D* BookPlateauMap(const TH3D* src, const std::string& name, const std::string& ztitle)
{
    const TAxis* ay = src->GetYaxis();   // pair pT
    const TAxis* az = src->GetZaxis();   // pair eta
    const int nx = ay->GetNbins(), ny = az->GetNbins();
    std::vector<double> xe(nx + 1), ye(ny + 1);
    for (int i = 0; i < nx; ++i) xe[i] = ay->GetBinLowEdge(i + 1);
    xe[nx] = ay->GetBinUpEdge(nx);
    for (int i = 0; i < ny; ++i) ye[i] = az->GetBinLowEdge(i + 1);
    ye[ny] = az->GetBinUpEdge(ny);
    auto* h = new TH2D(name.c_str(),
                       (";p_{T}^{pair} [GeV];#eta^{pair};" + ztitle).c_str(),
                       nx, xe.data(), ny, ye.data());
    h->SetDirectory(nullptr);
    return h;
}

// P is indexed [iy-1][iz-1] exactly as the tables build it (iy = pair-pT bin, iz = pair-eta bin).
// `sign` selects the series: "" = sign-integrated (nominal), "ss" = same sign, "os" = opposite
// sign. It only tags the key names, so one file holds all three and the fit stage picks a series
// by name -- no second file, no chance of pairing one sample's signs with another's.
void WritePlateauRootFile(const std::string& path, bool recreate, int step, const TH3D* src,
                          const std::vector<std::vector<PlateauCell>>& P,
                          double incl_mean, double incl_err, double incl_syst,
                          const std::string& sample, const std::string& wp_text,
                          const std::string& quantity, const std::string& sign = "")
{
    const std::string tag = "step" + std::to_string(step) + (sign.empty() ? "" : "_" + sign);
    TDirectory* prev = gDirectory;   // restored below: the caller keeps reading its input files
    TFile* f = TFile::Open(path.c_str(), recreate ? "RECREATE" : "UPDATE");
    if (!f || f->IsZombie())
        throw std::runtime_error("WritePlateauRootFile: cannot open " + path + " for writing");
    f->cd();

    TH2D* hval = BookPlateauMap(src, "h_" + tag + "_plateau",       quantity + " plateau");
    TH2D* herr = BookPlateauMap(src, "h_" + tag + "_plateau_err",   "stat. error on the plateau");
    TH2D* hrms = BookPlateauMap(src, "h_" + tag + "_plateau_rms",   "weighted RMS scatter");
    TH2D* hnb  = BookPlateauMap(src, "h_" + tag + "_plateau_nbins", "n #DeltaR bins in the window");
    TH2D* hsys = BookPlateauMap(src, "h_" + tag + "_plateau_syst",
                                Form("plateau-window systematic |p_{[%g,%g]} - p_{[%g,%g]}|",
                                     kPlateauLo, kPlateauHi, kPlateauSystLo, kPlateauSystHi));
    const int npt = (int)P.size(), neta = npt ? (int)P[0].size() : 0;
    for (int iy = 1; iy <= npt; ++iy) {
        for (int iz = 1; iz <= neta; ++iz) {
            const PlateauCell& p = P[iy - 1][iz - 1];
            // An unmeasurable cell (no usable dR bin in the window) is written as 0 with 0
            // error; the consumer MUST treat plateau <= 0 as "no plateau", never as a divisor.
            hval->SetBinContent(iy, iz, p.nb > 0 ? p.mean : 0.);
            hval->SetBinError  (iy, iz, p.nb > 0 ? p.err  : 0.);
            herr->SetBinContent(iy, iz, p.nb > 0 ? p.err  : 0.);
            hrms->SetBinContent(iy, iz, p.nb > 0 ? p.rms  : 0.);
            hnb ->SetBinContent(iy, iz, p.nb);
            // -1 = not evaluable (the alternative window had no usable bin), distinct from 0.
            hsys->SetBinContent(iy, iz, p.nb > 0 ? p.syst : -1.);
        }
    }
    // Inclusive (all cells) plateau, as a 1-bin histogram: value +- stat error.
    auto* hincl = new TH1D(("h_" + tag + "_plateau_inclusive").c_str(),
                           (";;" + quantity + " plateau (inclusive)").c_str(), 1, 0., 1.);
    hincl->SetDirectory(nullptr);
    hincl->SetBinContent(1, incl_mean);
    hincl->SetBinError(1, incl_err);
    auto* hisys = new TH1D(("h_" + tag + "_plateau_syst_inclusive").c_str(),
                           (";;" + quantity + " plateau-window systematic (inclusive)").c_str(),
                           1, 0., 1.);
    hisys->SetDirectory(nullptr);
    hisys->SetBinContent(1, incl_syst);

    for (TH1* h : {(TH1*)hval, (TH1*)herr, (TH1*)hrms, (TH1*)hnb, (TH1*)hsys,
                   (TH1*)hincl, (TH1*)hisys}) h->Write();

    // Provenance, so a stale file can be recognised without guessing.
    TNamed(("prov_" + tag).c_str(),
           Form("sample=%s; WP=%s; quantity=%s; plateau window dR in [%.2f,%.2f]; "
                "systematic window dR in [%.2f,%.2f]; source=plot_mc_trig_eff.cxx",
                sample.c_str(), wp_text.c_str(), quantity.c_str(),
                kPlateauLo, kPlateauHi, kPlateauSystLo, kPlateauSystHi)).Write();
    f->Close();
    if (prev) prev->cd();
    std::cout << "  wrote plateau map " << tag << " -> " << path << std::endl;
    delete hval; delete herr; delete hrms; delete hnb; delete hsys;
    delete hincl; delete hisys;
}

} // namespace

// ================================================================= main

void plot_mc_trig_eff(const std::string& sample = "pp", bool use_tight_wp = true)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gErrorIgnoreLevel = kWarning;

    const SampleCfg cfg = MakeCfg(sample, use_tight_wp);
    // Shared sample identity -- used here only to name the plateau ROOT file that the fit stage
    // reads back (dr_correction_sample_cfg.h).
    const DrCorrSample id = GetDrCorrSample(sample, use_tight_wp);

    // WP config (registry: Analysis/docs/muon_wp_registry.md): Tight nominal unsuffixed;
    // Medium inputs carry _medium_wp. BOTH the MC inputs AND the data tag-and-probe file are
    // WP-keyed (see MakeCfg) -- the data Medium reference is produced by running the data RDF
    // with isTight=false. Never compare across working points.
    const std::string wp_suf  = use_tight_wp ? "" : "_medium_wp";
    const std::string wp_text = use_tight_wp ? "Tight muons" : "Medium muons";
    const std::string headline = cfg.sample_text + ", " + wp_text;

    TFile* fmc   = OpenFile(cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf + ".root");
    TFile* fmc3  = OpenFile(cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf
                          + MCTrigEffPairPt::FileSuffix() + "_step3.root");
    TFile* ffit  = OpenFile(cfg.mc_dir + "single_mu_effcy_pT_fit_mc" + wp_suf + ".root");
    TFile* fdata = OpenFile(cfg.data_file);

    // Medium-WP plots go into a medium/ SUBDIRECTORY of each step dir (same filenames as
    // Tight); Tight (nominal) stays in the step dir root (user request, 2026-07-16).
    // The working point is in the TOP-LEVEL tree now (out_base = mc_based{_pt4bin}{_medium},
    // dr_correction_sample_cfg.h::DrCorrOutTag), so there is no per-step medium/ subdirectory.
    // This local copy kept creating one, leaving mc_based_medium/step3_dr_correction/medium/ --
    // the half-migrated layout the user asked to be rid of.
    const std::string wp_dir = "";
    const std::string dir1 = cfg.out_base + "step1_singles_data_mc/"  + wp_dir;
    // Step-2 is produced 3x (round-5 #3): the L1, HLT|L1 and full-chain efficiencies each go
    // into their OWN subdirectory of step2_dr_binned_singles/ (built inside the stage loop).
    const std::string dir2_base = cfg.out_base + "step2_dr_binned_singles/";
    const std::string dir3 = cfg.out_base + "step3_dr_correction/" + wp_dir;
    for (const auto& d : {dir1, dir3}) gSystem->mkdir(d.c_str(), kTRUE);

    // The two series are DIFFERENT estimators of the same efficiency, so each legend entry
    // states its own conditional probability in full (user decision 2026-08-04). No
    // abbreviations the audience has to decode ("T&P") and no drawing asides ("+ fit").
    const std::string mc_leg   = "MC, P(mu4 | reconstructed #mu)";
    const std::string data_leg = cfg.data_text +
                                 ", tag-and-probe P(2mu4 | mu4 tag, #DeltaR > 0.8)";

    // ================================================================
    // Step 1 (§3.1): data vs MC singles efficiency
    // ================================================================
    std::cout << "\n===== Step 1 (" << sample << ", " << wp_text << ") =====\n";

    struct Var { std::string mc, data, xtitle; double xlo, xhi; bool logx; };
    const std::vector<Var> vars = {
        {"pt",  "pt2nd",  "p_{T} [GeV]", 4.0, 60.0, true},
        {"eta", "eta2nd", "#eta",       -2.4,  2.4, false},
        {"phi", "phi2nd", "#phi",       -M_PI, M_PI, false}};

    for (const auto& v : vars) {
        TCanvas c(("c_step1_" + v.mc).c_str(), "", 1400, 800);
        c.Divide(2, 1);
        for (int ic = 0; ic < 2; ++ic) {
            c.cd(ic + 1);
            auto pads = SplitPadForRatio("s1_" + v.mc + "_" + std::to_string(ic), v.logx);

            TH1D* mnum = GetObj<TH1D>(fmc, "h_mc_" + v.mc + "_num_"   + kCharges[ic]);
            TH1D* mden = GetObj<TH1D>(fmc, "h_mc_" + v.mc + "_denom_" + kCharges[ic]);
            TH1D* dnum = GetObj<TH1D>(fdata, "h_" + v.data + cfg.ctr + "_" + kDataSigns[ic] + "_2mu4_sepr");
            TH1D* dden = GetObj<TH1D>(fdata, "h_" + v.data + cfg.ctr + "_" + kDataSigns[ic] + "_mu4_sepr");

            auto* gmc = BayesEff(mnum, mden);
            auto* gda = BayesEff(dnum, dden);
            StyleGraph(gmc, kMCColor, 21);
            StyleGraph(gda, kDataColor, 20);

            pads.first->cd();
            DrawEffFrame(v.xlo, v.xhi, "", 0.0, 1.1, "efficiency", true);
            DrawUnityLine(v.xlo, v.xhi);
            gda->Draw("PZ same");
            gmc->Draw("PZ same");

            DrawHeadline(headline + ", " + kChargeTex[ic], 0.14, 0.955, 0.05);
            // wide box + smaller text: the overlay data label ("... T&P P(2mu4 | mu4 tag,
            // DR>0.8)") is long and must not clip at the pad edge (plot review iter 1)
            auto* leg = new TLegend(0.18, 0.10, 0.93, 0.30);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.036);
            leg->AddEntry(gmc, mc_leg.c_str(), "lp");
            leg->AddEntry(gda, data_leg.c_str(), "lp");
            leg->Draw();

            // R3 ratio pad: data is the reference, so MC/data is the quantity of interest
            // (a flat offset = the known L1 over-efficiency; structure = something else).
            pads.second->cd();
            DrawRatioFrame(v.xlo, v.xhi, v.xtitle, "MC / data", 0.5, 2.2);
            auto* grat = DivideGraphClean(gmc, gda);
            StyleGraph(grat, kMCColor, 21, 0.9);
            grat->Draw("PZ same");
            MarkOffScale(grat, 0.5, 2.2, kMCColor);
            c.cd(ic + 1);

            // record MC/data ratio magnitudes (turn-on + plateau) from the pt hists
            if (v.mc == "pt") {
                std::cout << "  " << kCharges[ic] << " MC/data eff ratio vs pT:\n";
                for (double pt : {4.3, 5.5, 7.0, 10.0, 20.0, 40.0}) {
                    const int bm = mden->FindBin(pt);
                    const int bd = dden->FindBin(pt);
                    const double em = mden->GetBinContent(bm) > 0
                        ? mnum->GetBinContent(bm) / mden->GetBinContent(bm) : 0.;
                    const double ed = dden->GetBinContent(bd) > 0
                        ? dnum->GetBinContent(bd) / dden->GetBinContent(bd) : 0.;
                    printf("    pT=%5.1f  MC=%.4f  data=%.4f  MC/data=%.3f\n",
                           pt, em, ed, ed > 0 ? em / ed : 0.);
                }
            }
        }
        SaveCanvas(c, dir1 + "step1_eff_" + v.mc + ".png");
    }

    // --- pT in q.eta bins: per charge, 3x4 grid (10 bins + legend pad) ---
    // Optional comparison-MC file (cfg.cmp_fit_file, "noovl" only) -- opened once.
    TFile* fcmp = cfg.cmp_fit_file.empty() ? nullptr : OpenFile(cfg.cmp_fit_file);
    // Optional black-series-as-MC file (cfg.qeta_black_mc_file, "noovl" only): replaces the
    // data T&P black series on the q.eta panels with pp24-conditions MC.
    TFile* fblackmc = cfg.qeta_black_mc_file.empty() ? nullptr : OpenFile(cfg.qeta_black_mc_file);
    const bool qeta_all_mc = (fblackmc != nullptr);
    for (int ic = 0; ic < 2; ++ic) {
        TCanvas c(("c_step1_qeta_" + kCharges[ic]).c_str(), "", 1500, 2100);
        c.Divide(3, 4);  // ncols=3, nrows=4 (nrows >= ncols)
        for (size_t iq = 0; iq < kQEtaSuffix.size(); ++iq) {
            c.cd(static_cast<int>(iq) + 1);
            auto pads = SplitPadForRatio("s1q_" + kCharges[ic] + "_" + std::to_string(iq), true);

            auto* gmc = GetObj<TGraphAsymmErrors>(ffit,
                "g_mc_pt_vs_q_eta_" + kCharges[ic] + "_" + kQEtaSuffix[iq]);
            auto* fmy = GetObj<TF1>(ffit,
                "f_mc_pt_vs_q_eta_" + kCharges[ic] + "_" + kQEtaSuffix[iq]);
            // Black series: pp24-conditions MC (qeta_all_mc) OR data tag-and-probe.
            auto* gda = qeta_all_mc
                ? (TGraphAsymmErrors*)GetObj<TGraphAsymmErrors>(fblackmc,
                      "g_mc_pt_vs_q_eta_" + kCharges[ic] + "_" + kQEtaSuffix[iq])
                : GetObj<TGraphAsymmErrors>(fdata,
                      "g_pt2nd_vs_q_eta2nd" + cfg.ctr + "_" + kDataSigns[ic] +
                      "_2mu4_sepr_py_" + kQEtaSuffix[iq] + "_divided");

            gmc = (TGraphAsymmErrors*)gmc->Clone();
            gda = (TGraphAsymmErrors*)gda->Clone();
            // Strip the black graph's stored fit when it is MC (same TF1-auto-draw trap as the
            // comparison series): the black MC is shown as POINTS ONLY, no curve.
            if (qeta_all_mc) gda->GetListOfFunctions()->Clear();
            StyleGraph(gmc, kMCColor, 21, 0.7);
            StyleGraph(gda, kDataColor, 20, 0.7);

            pads.first->cd();
            gPad->SetLogx();
            DrawEffFrame(4.0, 60.0, "", 0.0, 1.1, "efficiency", true);
            DrawUnityLine(4.0, 60.0);
            auto* fdraw = (TF1*)fmy->Clone();
            fdraw->SetLineColor(kMCColor);
            fdraw->SetLineWidth(1);
            fdraw->Draw("same");
            gda->Draw("PZ same");
            // Comparison MC (no fit) UNDER the nominal MC, so the sample this plot set
            // belongs to stays on top and unobscured.
            if (fcmp) {
                auto* gcmp = (TGraphAsymmErrors*)GetObj<TGraphAsymmErrors>(fcmp,
                    "g_mc_pt_vs_q_eta_" + kCharges[ic] + "_" + kQEtaSuffix[iq])->Clone();
                // MUST strip the stored fit: FitMCSinglesEffcy writes these graphs AFTER
                // TGraph::Fit, so each carries its TF1 in its function list and ROOT draws
                // it automatically with "PZ same" -- silently adding a SECOND red curve
                // (FitMCSinglesEffcy sets the fit line red) that is neither requested nor
                // this sample's fit. The comparison sample is shown as POINTS ONLY.
                gcmp->GetListOfFunctions()->Clear();
                StyleGraph(gcmp, kCmpColor, 22, 0.7);
                gcmp->Draw("PZ same");
            }
            gmc->Draw("PZ same");

            TLatex tl;
            tl.SetNDC();
            tl.SetTextSize(0.075);
            tl.SetTextFont(42);
            tl.DrawLatex(0.30, 0.12, Form("%.1f < q#upoint#eta < %.1f",
                                          kQEtaRange[iq].first, kQEtaRange[iq].second));

            // R3 ratio pad: this-sample-MC / (black series) per q.eta bin. This is where the
            // (-2.4,-2.0) low-pT discrepancy lives -- the ratio makes it quantitative. When
            // the black series is pp24 MC (qeta_all_mc) the ratio IS the R10 MC-vs-MC
            // comparison (r17663 / pp24 cond.); otherwise it is the usual MC/data.
            pads.second->cd();
            gPad->SetLogx();
            DrawRatioFrame(4.0, 60.0, "p_{T} [GeV]",
                           qeta_all_mc ? "no overlay / pp cond." : "MC / data", 0.5, 2.6);
            auto* grat = DivideGraphClean(gmc, gda);
            StyleGraph(grat, kMCColor, 21, 0.7);
            grat->Draw("PZ same");
            MarkOffScale(grat, 0.5, 2.6, kMCColor);
            c.cd(static_cast<int>(iq) + 1);
        }
        // legend / label pad
        // legend/label pad = the first pad AFTER the q.eta panels. Derived, never hardcoded:
        // the bin count changed 11 -> 10 in round 8 and a literal 12 would have drawn the
        // legend into an occupied panel.
        c.cd(static_cast<int>(kQEtaSuffix.size()) + 1);
        DrawHeadline(headline + ", " + kChargeTex[ic], 0.02, 0.88, 0.048);
        auto* gm = new TGraphAsymmErrors(); StyleGraph(gm, kMCColor, 21);
        auto* gd = new TGraphAsymmErrors(); StyleGraph(gd, kDataColor, 20);
        if (qeta_all_mc) {
            // Three-way MC comparison (noovl): the legend spells out the reco-tag CONFIGURATION
            // that distinguishes the three curves, per user request. Smaller text so the full
            // labels are not clipped at the legend-pad width.
            // Colour -> reco-tag CONFIGURATION (user's exact wording). r-tag names are dropped
            // from the legend text (they clip at the pad width) and given in the note below;
            // SetMargin trims the symbol column so the full condition strings fit.
            auto* leg = new TLegend(0.02, 0.50, 0.98, 0.74);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.030);
            leg->SetMargin(0.12);
            leg->AddEntry(gm, "MC, Pb+Pb 2023 conditions, no overlay", "lp");
            leg->AddEntry(gd, "MC, pp 2024 conditions", "lp");
            if (fcmp) {
                auto* gc = new TGraphAsymmErrors(); StyleGraph(gc, kCmpColor, 22);
                leg->AddEntry(gc, "MC, Pb+Pb 2023 conditions, with HIJING overlay", "lp");
            }
            leg->Draw();
            // No explanatory prose on the canvas: the three samples are identified by the
            // legend, and the ratio pad by its own y-axis title. What the samples mean and why
            // they are compared belongs in the tracking doc, not in front of the audience.
        } else {
            // Each entry carries its OWN definition, so no separate explanatory block is needed
            // (0.034, not 0.05: the full conditional probabilities must fit the pad width).
            auto* leg = new TLegend(0.02, 0.45, 0.98, 0.80);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            // 0.026 + a trimmed symbol column: the data entry is ~50 glyphs and was clipped at
            // the pad edge, cutting off the end of its DeltaR condition.
            leg->SetTextSize(0.026);
            leg->SetMargin(0.12);
            leg->AddEntry(gm, mc_leg.c_str(), "lp");
            leg->AddEntry(gd, data_leg.c_str(), "lp");
            if (fcmp) {
                auto* gc = new TGraphAsymmErrors(); StyleGraph(gc, kCmpColor, 22);
                leg->AddEntry(gc, cfg.cmp_text.c_str(), "lp");
            }
            leg->Draw();
        }
        SaveCanvas(c, dir1 + "step1_eff_pt_in_q_eta_bins_" + kCharges[ic] + ".png");
    }

    // ================================================================
    // Step 2 (§3.2): DeltaR-binned MC singles efficiency
    // ================================================================
    std::cout << "\n===== Step 2 (" << sample << ", " << wp_text << ") =====\n";

    const std::string leg_incl = "all #DeltaR";

    // --- 1D pt and 1D q.eta, 2-pad canvases (mu+ | mu-) -------------
    struct PairVar { std::string tag, xtitle; double xlo, xhi; bool logx; };
    const std::vector<PairVar> pvars = {
        {"pt",    "p_{T} [GeV]",   4.0, 60.0, true},
        {"q_eta", "q#upoint#eta", -2.4,  2.4, false}};

    // round-5 #3: produce Step-2 THREE times -- the L1 leg, the HLT|L1 leg, and the full mu4
    // chain -- each into its own subdirectory. (numk/denk) pick the numerator/denominator hist
    // infix: full_chain = num/denom, L1 = numl1/denom, HLT|L1 = numhlt/numl1, so
    // eff(chain) = eff(L1) * eff(HLT|L1). The L1/HLT hists are empty on pre-reskim NTUPs (no
    // per-muon L1 branch) -> those subdirs are only meaningful after the re-skim.
    // `sub` names the output SUBDIRECTORY only; what the audience sees is `eff`, the efficiency
    // itself, drawn as the legend header. The subdirectory token is never put on a canvas.
    struct Step2Stage { std::string sub, numk, denk, eff; };
    const std::vector<Step2Stage> step2_stages = {
        {"full_chain", "num",    "denom", "P(mu4 | reconstructed #mu)"},
        {"L1",         "numl1",  "denom", "P(L1 MU3V | reconstructed #mu)"},
        {"HLT",        "numhlt", "numl1", "P(mu4 HLT | L1 MU3V, reconstructed #mu)"}};
    for (const auto& st : step2_stages) {
      const std::string dir2 = dir2_base + st.sub + "/" + wp_dir;
      gSystem->mkdir(dir2.c_str(), kTRUE);
      const std::string NUMK = st.numk, DENK = st.denk;
      std::cout << "  -- Step-2 stage: " << st.sub << " (" << st.eff << ")\n";

    for (const auto& v : pvars) {
        TCanvas c(("c_step2_" + v.tag).c_str(), "", 1400, 800);
        c.Divide(2, 1);
        for (int ic = 0; ic < 2; ++ic) {
            c.cd(ic + 1);
            auto pads = SplitPadForRatio("s2_" + v.tag + "_" + std::to_string(ic), v.logx);
            pads.first->cd();

            DrawEffFrame(v.xlo, v.xhi, "", 0.0, 1.1, "efficiency", true);
            DrawUnityLine(v.xlo, v.xhi);

            // coarse edges for this variable (empty => keep the native binning)
            const std::vector<double> cedges =
                !cfg.step2_coarse ? std::vector<double>{}
                : (v.tag == "pt" ? kStep2CoarsePtTarget : Step2CoarseQEtaTarget());
            auto Coarsen = [&](TH1* h, const std::string& nm) -> TH1* {
                return cedges.empty() ? h : RebinTo(h, cedges, nm);
            };

            // Step-1 inclusive singles reference (thin black). No 1D q.eta singles hist
            // exists -> project the singles 2D (x = q.eta) over all pT.
            TH1* rnum = nullptr;
            TH1* rden = nullptr;
            if (v.tag == "pt") {
                rnum = GetObj<TH1D>(fmc, "h_mc_pt_" + NUMK + "_" + kCharges[ic]);
                rden = GetObj<TH1D>(fmc, "h_mc_pt_" + DENK + "_" + kCharges[ic]);
            } else {
                TH2D* h2n = GetObj<TH2D>(fmc, "h_mc_pt_vs_q_eta_" + NUMK + "_" + kCharges[ic]);
                TH2D* h2d = GetObj<TH2D>(fmc, "h_mc_pt_vs_q_eta_" + DENK + "_" + kCharges[ic]);
                rnum = h2n->ProjectionX(Form("px_num_%s_%d", v.tag.c_str(), ic));
                rden = h2d->ProjectionX(Form("px_den_%s_%d", v.tag.c_str(), ic));
            }
            rnum = Coarsen(rnum, Form("rb_ref_n_%s_%d", v.tag.c_str(), ic));
            rden = Coarsen(rden, Form("rb_ref_d_%s_%d", v.tag.c_str(), ic));
            auto* gref = BayesEff(rnum, rden);
            StyleGraph(gref, kBlack, 1, 0.4, 1);

            auto* leg = new TLegend(0.38, 0.10, 0.92, 0.36);
            leg->SetBorderSize(0);
            // Semi-opaque backing so the legend stays readable inside a dense error-bar
            // cloud -- but it is drawn FIRST, BEFORE the data, so the fill can never wash
            // out real points (it did, in the endcap of the q.eta panels: review iter 1).
            leg->SetFillColor(kWhite);   // SOLID: alpha does not render in batch PNG
            leg->SetFillStyle(1001);
            leg->SetTextSize(0.042);

            std::vector<TGraphAsymmErrors*> gdr;
            for (size_t id = 0; id < kDrSuffix.size(); ++id) {
                TH1* n = GetObj<TH1D>(fmc, "h_mc_pair_" + v.tag + "_" + NUMK + "_" +
                                            kCharges[ic] + "_" + kDrSuffix[id]);
                TH1* d = GetObj<TH1D>(fmc, "h_mc_pair_" + v.tag + "_" + DENK + "_" +
                                            kCharges[ic] + "_" + kDrSuffix[id]);
                n = Coarsen(n, Form("rb_n_%s_%d_%zu", v.tag.c_str(), ic, id));
                d = Coarsen(d, Form("rb_d_%s_%d_%zu", v.tag.c_str(), ic, id));
                auto* g = BayesEff(n, d);
                StyleGraph(g, kDrColor[id], kDrMarker[id], 0.8);
                leg->AddEntry(g, kDrTex[id].c_str(), "lp");
                gdr.push_back(g);
            }
            leg->AddEntry(gref, leg_incl.c_str(), "l");
            leg->SetHeader(st.eff.c_str());    // which efficiency this stage plots (round-5 #3)
            leg->Draw();                       // legend first ...
            gref->Draw("LX same");             // ... then the data on top of it
            for (auto* g : gdr) g->Draw("PZ same");
            DrawHeadline(headline + ", " + kChargeTex[ic], 0.14, 0.955, 0.05);

            // R3 ratio pad: §3.2 asks whether the three DeltaR series AGREE at fixed
            // kinematics. Ratio to the isolated (DeltaR >= 1.0) series = that test.
            pads.second->cd();
            DrawRatioFrame(v.xlo, v.xhi, v.xtitle, "ratio to #DeltaR #geq 1", 0.4, 1.9);
            for (size_t id = 0; id + 1 < gdr.size(); ++id) {
                auto* g = DivideGraphClean(gdr[id], gdr.back());
                StyleGraph(g, kDrColor[id], kDrMarker[id], 0.8);
                g->Draw("PZ same");
                MarkOffScale(g, 0.4, 1.9, kDrColor[id]);
            }
            c.cd(ic + 1);
        }
        SaveCanvas(c, dir2 + "step2_eff_" + v.tag + "_dr_bins.png");
    }

    // --- pT in q.eta bins per charge, 3 DeltaR lines per pad --------
    for (int ic = 0; ic < 2; ++ic) {
        TCanvas c(("c_step2_qeta_" + kCharges[ic]).c_str(), "", 1500, 2100);
        c.Divide(3, 4);
        for (size_t iq = 0; iq < kQEtaSuffix.size(); ++iq) {
            c.cd(static_cast<int>(iq) + 1);
            auto pads = SplitPadForRatio("s2q_" + kCharges[ic] + "_" + std::to_string(iq), true);
            pads.first->cd();
            gPad->SetLogx();
            DrawEffFrame(4.0, 60.0, "", 0.0, 1.1, "efficiency", true);
            DrawUnityLine(4.0, 60.0);
            std::vector<TGraphAsymmErrors*> gdr;
            for (size_t id = 0; id < kDrSuffix.size(); ++id) {
                TH2D* h2n = GetObj<TH2D>(fmc, "h_mc_pair_pt_vs_q_eta_" + NUMK + "_" +
                                              kCharges[ic] + "_" + kDrSuffix[id]);
                TH2D* h2d = GetObj<TH2D>(fmc, "h_mc_pair_pt_vs_q_eta_" + DENK + "_" +
                                              kCharges[ic] + "_" + kDrSuffix[id]);
                const int blo = h2n->GetXaxis()->FindBin(kQEtaRange[iq].first  + 1e-6);
                const int bhi = h2n->GetXaxis()->FindBin(kQEtaRange[iq].second - 1e-6);
                TH1* n = h2n->ProjectionY(Form("py_n_%d_%zu_%zu", ic, iq, id), blo, bhi);
                TH1* d = h2d->ProjectionY(Form("py_d_%d_%zu_%zu", ic, iq, id), blo, bhi);
                if (cfg.step2_coarse) {
                    n = RebinTo(n, kStep2CoarsePtTarget, Form("rb_qn_%d_%zu_%zu", ic, iq, id));
                    d = RebinTo(d, kStep2CoarsePtTarget, Form("rb_qd_%d_%zu_%zu", ic, iq, id));
                }
                auto* g = BayesEff(n, d);
                StyleGraph(g, kDrColor[id], kDrMarker[id], 0.7);
                g->Draw("PZ same");
                gdr.push_back(g);
            }
            TLatex tl;
            tl.SetNDC();
            tl.SetTextSize(0.075);
            tl.SetTextFont(42);
            tl.DrawLatex(0.30, 0.12, Form("%.1f < q#upoint#eta < %.1f",
                                          kQEtaRange[iq].first, kQEtaRange[iq].second));

            // R3 ratio pad: the §3.2 factorization test = do the close-pair series agree
            // with the isolated (DeltaR >= 1.0) one at fixed kinematics?
            pads.second->cd();
            gPad->SetLogx();
            DrawRatioFrame(4.0, 60.0, "p_{T} [GeV]", "ratio to #DeltaR #geq 1", 0.4, 1.9);
            for (size_t id = 0; id + 1 < gdr.size(); ++id) {
                auto* g = DivideGraphClean(gdr[id], gdr.back());
                StyleGraph(g, kDrColor[id], kDrMarker[id], 0.7);
                g->Draw("PZ same");
                MarkOffScale(g, 0.4, 1.9, kDrColor[id]);
            }
            c.cd(static_cast<int>(iq) + 1);
        }
        // legend/label pad = the first pad AFTER the q.eta panels. Derived, never hardcoded:
        // the bin count changed 11 -> 10 in round 8 and a literal 12 would have drawn the
        // legend into an occupied panel.
        c.cd(static_cast<int>(kQEtaSuffix.size()) + 1);
        DrawHeadline(headline + ", " + kChargeTex[ic], 0.02, 0.88, 0.048); // 0.048: long overlay headline + charge must fit (review iter 1)
        auto* leg = new TLegend(0.05, 0.35, 0.95, 0.80);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.06);
        leg->SetHeader(st.eff.c_str());   // which efficiency this stage plots (round-5 #3)
        for (size_t id = 0; id < kDrSuffix.size(); ++id) {
            auto* g = new TGraphAsymmErrors();
            StyleGraph(g, kDrColor[id], kDrMarker[id]);
            leg->AddEntry(g, kDrTex[id].c_str(), "lp");
        }
        leg->Draw();
        SaveCanvas(c, dir2 + "step2_eff_pt_in_q_eta_bins_" + kCharges[ic] +
                          "_dr_bins.png");
    }
    }  // end round-5 #3 Step-2 stage loop (full_chain / L1 / HLT)

    // ================================================================
    // Step-1 SANITY CHECK (§3.5, round 7): is the MC >> data forward single-muon efficiency
    // caused by "bad" muons/events rather than by the simulation?
    //   variant orig = the round-7 nominal selection
    //           vtx  = + exactly ONE reconstructed track-bearing primary vertex (no pile-up)
    //           ptm  = + |truth pT - reco pT| / truth pT < threshold (no badly measured muons)
    //           both = + both
    // The overlay is compared to the ORIGINAL MC, never to data: neither extra requirement has
    // a data analogue (data has no truth, and the vertex requirement would change the event
    // sample rather than the muon selection). Graceful skip if the _sanity.root is absent.
    // ================================================================
    {
        const std::string sanity_path =
            cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf + "_sanity.root";
        TFile* fsan = TFile::Open(sanity_path.c_str(), "READ");
        if (!fsan || fsan->IsZombie()) {
            std::cout << "\n[Step-1 sanity] no " << sanity_path
                      << " -- skipping (run FillMCTrigEffHists do_sanity=true first).\n";
        } else {
            std::cout << "\n===== Step-1 SANITY CHECK (" << sample << ", " << wp_text << ") =====\n";
            const std::string dirs = cfg.out_base + "step1_sanity_check/" + wp_dir;
            gSystem->mkdir(dirs.c_str(), kTRUE);

            // `key` is the internal variant token (histogram names only); `tex` is what the
            // audience reads and must state the requirement IN FULL, including the numerical
            // threshold -- taken from the shared header, never retyped.
            struct Var { std::string key, tex; Color_t col; Style_t mk; };
            std::vector<Var> vars = {
                {"orig", "nominal muon selection",                          kBlack,     20},
                {"vtx",  "+ exactly 1 reconstructed vertex",                kBlue + 1,  21},
                {"ptm",  Form("+ |#Deltap_{T}|/p_{T}^{truth} < %.2f",
                              kSanityPtMatchThr),                           kGreen + 2, 22},
                {"both", "+ both requirements",                             kRed + 1,   23}};

            // A variant with an EMPTY denominator is not a measurement of zero -- it means the
            // requirement selected nothing in this sample. Drop those series and say so on the
            // canvas instead of drawing empty graphs. (Note the DIFFERENT, benign case: in the
            // HIJING overlay every event has exactly one track-bearing vertex, so the `vtx`
            // variant is identical to `orig` rather than empty -- a no-op, not a failure.)
            std::vector<std::string> dropped;
            {
                std::vector<Var> keep;
                for (const auto& v : vars) {
                    double d = 0;
                    for (const auto& chg : kCharges)
                        d += GetObj<TH1D>(fsan, "h_sanity_pt_denom_" + chg + "_" + v.key)->Integral();
                    if (d > 0) keep.push_back(v);
                    else { // record the READABLE requirement, not the internal variant token:
                           // "vtx"/"ptm" mean nothing to the audience. Strip the leading "+ ".
                           dropped.push_back(v.tex.rfind("+ ", 0) == 0 ? v.tex.substr(2) : v.tex);
                           std::cout << "  [sanity] variant '" << v.key
                                     << "' has an EMPTY denominator -> requirement inapplicable "
                                        "to this sample; series dropped\n"; }
                }
                vars = keep;
            }
            const std::string drop_note = dropped.empty() ? std::string()
                : ("no muon in this sample satisfies: " +
                   [&]{ std::string s; for (size_t i = 0; i < dropped.size(); ++i)
                        s += (i ? "; " : "") + dropped[i]; return s; }());

            auto eff_of = [&](const std::string& base, const std::string& chg,
                              const std::string& v) -> TGraphAsymmErrors* {
                TH1D* n = GetObj<TH1D>(fsan, "h_sanity_" + base + "_num_"   + chg + "_" + v);
                TH1D* d = GetObj<TH1D>(fsan, "h_sanity_" + base + "_denom_" + chg + "_" + v);
                return BayesEff(n, d);
            };

            // ---- (a) eff vs pT and vs q.eta, one pad per charge, ratio pad vs `orig` ----
            struct Obs { std::string base, xt; double xlo, xhi; bool logx; };
            const std::vector<Obs> obs = {{"pt", "p_{T} [GeV]", 4.0, 60.0, true},
                                          {"q_eta", "q#upoint#eta", -2.4, 2.4, false}};
            for (const auto& O : obs) {
                TCanvas c(("c_sanity_" + O.base).c_str(), "", 1500, 700);
                c.Divide(2, 1);
                for (int ic = 0; ic < 2; ++ic) {
                    c.cd(ic + 1);
                    auto pads = SplitPadForRatio("san_" + O.base + "_" + kCharges[ic], O.logx);
                    std::vector<TGraphAsymmErrors*> gs;
                    for (const auto& v : vars) {
                        auto* g = eff_of(O.base, kCharges[ic], v.key);
                        StyleGraph(g, v.col, v.mk, 0.8);
                        gs.push_back(g);
                    }
                    pads.first->cd();
                    if (O.logx) gPad->SetLogx();
                    DrawEffFrame(O.xlo, O.xhi, "", 0.0, 1.1, "efficiency", true);
                    DrawUnityLine(O.xlo, O.xhi);
                    for (auto* g : gs) g->Draw("PZ same");
                    auto* leg = new TLegend(0.40, 0.15, 0.93, 0.42);
                    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
                    for (size_t i = 0; i < vars.size(); ++i)
                        leg->AddEntry(gs[i], vars[i].tex.c_str(), "lp");
                    leg->Draw();
                    TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.045);
                    tl.DrawLatex(0.18, 0.88, kChargeTex[ic].c_str());

                    pads.second->cd();
                    if (O.logx) gPad->SetLogx();
                    DrawEffFrame(O.xlo, O.xhi, O.xt, 0.90, 1.10, "requirement / nominal");
                    DrawUnityLine(O.xlo, O.xhi);
                    for (size_t i = 1; i < gs.size(); ++i) {   // skip `orig` (ratio 1 by construction)
                        auto* gr = DivideGraphClean(gs[i], gs[0]);
                        StyleGraph(gr, vars[i].col, vars[i].mk, 0.8);
                        gr->Draw("PZ same");
                    }
                }
                // c.cd(0) FIRST: the loop above leaves gPad pointing at the RATIO sub-pad of the
                // second charge, so the headline used to be drawn inside it -- rendered tiny and
                // overlapping the ratio axis labels on every sanity canvas.
                c.cd(0);
                DrawHeadline(headline);
                if (!drop_note.empty()) {
                    TLatex nt; nt.SetNDC(); nt.SetTextFont(42); nt.SetTextSize(0.022);
                    nt.SetTextColor(kGray + 3); nt.DrawLatex(0.06, 0.015, drop_note.c_str());
                }
                SaveCanvas(c, dirs + "sanity_eff_" + O.base + ".png");
            }

            // ---- (b) pT turn-on per fine q.eta bin (the forward bins are the whole point) ----
            for (int ic = 0; ic < 2; ++ic) {
                TCanvas c(("c_sanity_qeta_" + kCharges[ic]).c_str(), "", 1500, 2100);
                c.Divide(3, 4);
                for (size_t iq = 0; iq < kQEtaSuffix.size(); ++iq) {
                    c.cd(static_cast<int>(iq) + 1);
                    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetLogx();
                    DrawEffFrame(4.0, 60.0, "p_{T} [GeV]", 0.0, 1.15);
                    DrawUnityLine(4.0, 60.0);
                    std::vector<TGraphAsymmErrors*> gs;
                    for (const auto& v : vars) {
                        TH2D* h2n = GetObj<TH2D>(fsan, "h_sanity_pt_vs_q_eta_num_"   + kCharges[ic] + "_" + v.key);
                        TH2D* h2d = GetObj<TH2D>(fsan, "h_sanity_pt_vs_q_eta_denom_" + kCharges[ic] + "_" + v.key);
                        const int b1 = h2n->GetXaxis()->FindBin(kQEtaRange[iq].first  + 1e-6);
                        const int b2 = h2n->GetXaxis()->FindBin(kQEtaRange[iq].second - 1e-6);
                        TH1D* n = h2n->ProjectionY(Form("san_n_%d_%zu_%s", ic, iq, v.key.c_str()), b1, b2);
                        TH1D* d = h2d->ProjectionY(Form("san_d_%d_%zu_%s", ic, iq, v.key.c_str()), b1, b2);
                        auto* g = BayesEff(n, d);
                        StyleGraph(g, v.col, v.mk, 0.7);
                        g->Draw("PZ same");
                        gs.push_back(g);
                        delete n; delete d;
                    }
                    TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.050);
                    tl.DrawLatex(0.20, 0.90, Form("%.2f < q#upoint#eta < %.2f",
                                                  kQEtaRange[iq].first, kQEtaRange[iq].second));
                    if (iq == 0) {
                        auto* leg = new TLegend(0.30, 0.14, 0.95, 0.40);
                        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
                        for (size_t i = 0; i < vars.size(); ++i)
                            leg->AddEntry(gs[i], vars[i].tex.c_str(), "lp");
                        leg->Draw();
                    }
                }
                c.cd(0);                       // same sub-pad trap as the (a) canvases above
                // RAW TLatex, not DrawHeadline: on the CANVAS the text size is a fraction of the
                // 2100 px canvas height, and DrawHeadline's 0.030 floor (meant for a pad) would
                // render an 60 px headline that runs off the right edge.
                TLatex hl; hl.SetNDC(); hl.SetTextFont(42); hl.SetTextSize(0.012);
                hl.DrawLatex(0.02, 0.988, (headline + ", " + kChargeTex[ic]).c_str());
                if (!drop_note.empty()) {
                    TLatex nt; nt.SetNDC(); nt.SetTextFont(42); nt.SetTextSize(0.014);
                    nt.SetTextColor(kGray + 3); nt.DrawLatex(0.04, 0.008, drop_note.c_str());
                }
                SaveCanvas(c, dirs + "sanity_eff_pt_in_q_eta_bins_" + kCharges[ic] + ".png");
            }

            // ---- (c) pass fractions + the forward-bin verdict numbers, as a table ----
            {
                std::ofstream os(dirs + "sanity_pass_fractions.txt");
                os << "# Step-1 sanity check (mc_trigger_efficiency.md §3.5)\n";
                os << "# sample=" << sample << "  WP=" << wp_text << "  " << cfg.sample_text << "\n";
                if (!drop_note.empty())
                    os << "# NOTE: " << drop_note
                       << "  (an empty denominator means the requirement cannot be applied to this\n"
                          "#       sample -- e.g. the HIJING overlay reconstructs no track-bearing\n"
                          "#       primary vertex at all -- NOT that the efficiency is zero)\n";
                os << "# Fraction of the round-7-selected MC muons passing each extra requirement\n";
                os << "# (weighted by the MC event weight, and raw), plus the integrated eps(mu4).\n\n";
                const double wall = GetObj<TH1D>(fsan, "h_sanity_count_orig")->Integral();
                const double nall = GetObj<TH1D>(fsan, "h_sanity_rawcount_orig")->Integral();
                os << std::left << std::setw(34) << "requirement" << std::setw(14) << "weighted %"
                   << std::setw(14) << "raw %" << std::setw(16) << "raw N"
                   << std::setw(12) << "eps(mu+)" << std::setw(12) << "eps(mu-)" << "\n";
                for (const auto& v : vars) {
                    const double wv = GetObj<TH1D>(fsan, "h_sanity_count_"    + v.key)->Integral();
                    const double nv = GetObj<TH1D>(fsan, "h_sanity_rawcount_" + v.key)->Integral();
                    auto eff = [&](const std::string& chg) {
                        TH1D* n = GetObj<TH1D>(fsan, "h_sanity_pt_num_"   + chg + "_" + v.key);
                        TH1D* d = GetObj<TH1D>(fsan, "h_sanity_pt_denom_" + chg + "_" + v.key);
                        return d->Integral() > 0 ? n->Integral() / d->Integral() : -1.0;
                    };
                    os << std::left << std::setw(34) << v.key
                       << std::setw(14) << Form("%.4f", wall > 0 ? 100 * wv / wall : -1.)
                       << std::setw(14) << Form("%.4f", nall > 0 ? 100 * nv / nall : -1.)
                       << std::setw(16) << Form("%.0f", nv)
                       << std::setw(12) << Form("%.4f", eff("muplus"))
                       << std::setw(12) << Form("%.4f", eff("muminus")) << "\n";
                }
                os << "\n# eps(mu4) in the two split forward q.eta bins, pT 4-6 GeV "
                      "(the bins where MC >> data):\n";
                os << std::left << std::setw(28) << "q.eta bin";
                for (const auto& v : vars) os << std::setw(22) << v.key;
                os << "\n";
                for (int iq = 0; iq < 2; ++iq) {          // (-2.4,-2.2) and (-2.2,-2.0)
                    os << std::left << std::setw(28)
                       << Form("[%.1f,%.1f) mu+/mu-", kQEtaRange[iq].first, kQEtaRange[iq].second);
                    for (const auto& v : vars) {
                        double nn = 0, dd = 0;
                        for (int ic = 0; ic < 2; ++ic) {
                            TH2D* h2n = GetObj<TH2D>(fsan, "h_sanity_pt_vs_q_eta_num_"   + kCharges[ic] + "_" + v.key);
                            TH2D* h2d = GetObj<TH2D>(fsan, "h_sanity_pt_vs_q_eta_denom_" + kCharges[ic] + "_" + v.key);
                            const int b1 = h2n->GetXaxis()->FindBin(kQEtaRange[iq].first  + 1e-6);
                            const int b2 = h2n->GetXaxis()->FindBin(kQEtaRange[iq].second - 1e-6);
                            const int p1 = h2n->GetYaxis()->FindBin(4.0 + 1e-6);
                            const int p2 = h2n->GetYaxis()->FindBin(6.0 - 1e-6);
                            nn += h2n->Integral(b1, b2, p1, p2);
                            dd += h2d->Integral(b1, b2, p1, p2);
                        }
                        os << std::setw(22) << Form("%.4f", dd > 0 ? nn / dd : -1.);
                    }
                    os << "\n";
                }
                std::cout << "  wrote " << dirs << "sanity_pass_fractions.txt\n";
            }
            fsan->Close();
        }
    }

    // ================================================================
    // Step 3 (§3.3): eps_dR(dR) = inverse-weighted num / denom
    // ================================================================
    std::cout << "\n===== Step 3 (" << sample << ", " << wp_text << ") =====\n";

    // Central value = num/denom; ERROR = the conditional (binomial-correct) form, NOT
    // TH1::Divide's independent propagation (see SetConditionalRatioErrors).
    auto MakeRatio = [&](const std::string& tag) -> TH1D* {
        TH1D* num = GetObj<TH1D>(fmc3, "h_mc_dr_" + tag + "_num");
        TH1D* den = GetObj<TH1D>(fmc3, "h_mc_dr_" + tag + "_denom");
        TH1D* eA  = GetObj<TH1D>(fmc3, "h_mc_dr_" + tag + "_errA");
        TH1D* eB  = GetObj<TH1D>(fmc3, "h_mc_dr_" + tag + "_errB");
        auto* r = (TH1D*)num->Clone(("r_dr_" + tag).c_str());
        r->SetDirectory(nullptr);
        r->Divide(den);
        SetConditionalRatioErrors(r, den, eA, eB);
        return r;
    };
    TH1D* r_zoom = MakeRatio("zoom");
    TH1D* r_full = MakeRatio("full");

    const auto plateau = PlateauWeightedMean(r_full, kPlateauLo, kPlateauHi);
    const auto plateau_alt = PlateauWeightedMean(r_full, kPlateauSystLo, kPlateauSystHi);
    const double plateau_syst = std::fabs(plateau.first - plateau_alt.first);
    printf("  %s plateau (weighted mean, dR in [%g,%g]): %.4f +- %.4f"
           "   [window syst vs dR in [%g,%g] = %.4f]\n",
           sample.c_str(), kPlateauLo, kPlateauHi, plateau.first, plateau.second,
           kPlateauSystLo, kPlateauSystHi, plateau_syst);

    // REMOVED (round 9, user request): the two pair-pT/pair-eta INTEGRATED canvases
    // (step3_eps_dr_zoom.png, step3_eps_dr_full.png). An inclusive curve averages the correction
    // over the whole (pair pT, pair eta) plane, which is not a quantity the analysis applies --
    // the correction is used per cell. The inclusive plateau itself is still MEASURED above and
    // still written to the plateau ROOT file (h_step3_plateau_inclusive), which the fit stage
    // and the guard both read; only the two figures are gone. The per-cell distributions now
    // live under step3_dr_correction/<dR range>/, one PNG per pair-pT bin (see below).

    // --- pair-pT-binned zoom ratio, one series per COARSE pair-pT bin ---
    // Round 7: this panel used to slice a SEPARATE dR x fine-pair-pT 2D (pT_bins_120, grouped
    // 8-13.8/13.8-23.6/23.6-40.6/40.6-120) -- a second pair-pT binning inconsistent with the
    // canonical ParamsSet::pair_pt_coarse_bins used by the pair-eta panels, the
    // plateau tables, Step 4 and crossx. It now projects the SAME 3D as everything else,
    // integrating over ALL pair-eta, so 1D/2D/3D cannot disagree (see CLAUDE.md: pair-pT and
    // pair-eta binnings come from ParamsSet.h and are never re-invented per plot).
    {
        TH3D* h3zn_s = GetObj<TH3D>(fmc3, "h_mc_dr_zoom_vs_pt_eta_num");
        TH3D* h3zd_s = GetObj<TH3D>(fmc3, "h_mc_dr_zoom_vs_pt_eta_denom");
        TH3D* h3zA_s = GetObj<TH3D>(fmc3, "h_mc_dr_zoom_vs_pt_eta_errA");
        TH3D* h3zB_s = GetObj<TH3D>(fmc3, "h_mc_dr_zoom_vs_pt_eta_errB");
        const int npt_s  = h3zn_s->GetYaxis()->GetNbins();
        const int neta_s = h3zn_s->GetZaxis()->GetNbins();
        // last slice: bright kMagenta, NOT kMagenta+2 -- the darkened shade reads as another
        // dark red/blue against kRed+1 / kBlue+1 and the series cannot be told apart
        // 8 pair-pT bins since round 8: these are indexed by pair-pT bin BELOW WITHOUT a
        // modulo, so a 4-entry palette was an out-of-bounds read (undefined behaviour), not
        // just a colour clash. One distinct colour+marker per bin, and the indexing is
        // guarded with % anyway so a future bin-count change cannot resurrect the bug.
        const std::vector<Color_t> scol   = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta,
                                            kOrange + 7, kCyan + 2, kViolet + 1, kBlack};
        const std::vector<Style_t> smark  = {20, 21, 22, 23, 33, 34, 29, 24};

        std::vector<TH1D*> ratios;
        std::vector<std::string> slabels;
        for (int iy = 1; iy <= npt_s; ++iy) {
            TH1D* n = h3zn_s->ProjectionX(Form("s3_n_%d", iy), iy, iy, 1, neta_s, "e");
            TH1D* d = h3zd_s->ProjectionX(Form("s3_d_%d", iy), iy, iy, 1, neta_s, "e");
            TH1D* a = h3zA_s->ProjectionX(Form("s3_a_%d", iy), iy, iy, 1, neta_s, "e");
            TH1D* b = h3zB_s->ProjectionX(Form("s3_b_%d", iy), iy, iy, 1, neta_s, "e");
            auto* r = (TH1D*)n->Clone(Form("s3_r_%d", iy));
            r->SetDirectory(nullptr);
            r->Divide(d);
            SetConditionalRatioErrors(r, d, a, b);
            delete a; delete b;
            ratios.push_back(r);
            const double plo = h3zn_s->GetYaxis()->GetBinLowEdge(iy);
            const double phi = h3zn_s->GetYaxis()->GetBinUpEdge(iy);
            slabels.push_back(Form("%.1f < p_{T}^{pair} < %.1f GeV", plo, phi));
        }
        
        // ---- SPLIT INTO TWO PNGs (round 8, user) -------------------------------------------
        // With 8 pair-pT bins, overlaying every series in one pad is unreadable -- the curves sit
        // on top of each other and the small-DeltaR structure the panel exists to show is lost.
        // Draw the LOWER half of the pair-pT bins in one file and the UPPER half in another.
        // The y-range is computed PER HALF so each file is scaled to the data it actually shows.
        const int nhalf = (npt_s + 1) / 2;
        for (int half = 0; half < 2; ++half) {
            const int is_lo = half * nhalf;
            const int is_hi = std::min(npt_s, (half + 1) * nhalf);
            if (is_lo >= is_hi) continue;

            TCanvas c(Form("c_step3_ptslices_%d", half), "", 900, 700);
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);
            // Reserved strip ABOVE the frame for the legend: with the series filling the pad,
            // every in-frame position lands on data and a white legend backing would only hide
            // the points behind it.
            gPad->SetTopMargin(0.22);

            double ymax = 0.;
            for (int is = is_lo; is < is_hi; ++is)
                for (int i = 1; i <= ratios[is]->GetNbinsX(); ++i)
                    ymax = std::max(ymax, ratios[is]->GetBinContent(i) + ratios[is]->GetBinError(i));
            // CAP the auto-range: max+error is set by the noisiest bin, which would squeeze the
            // structure into the bottom sliver. Every point pushed off scale is MARKED with an
            // arrow and listed -- silently dropping data from a physics figure is not allowed.
            ymax = std::min(std::max(1.15, 1.15 * ymax), 3.0);
            DrawEffFrame(0.0, 1.0, "#DeltaR", 0.0, ymax, cfg.eps_dr_text);
            DrawUnityLine(0.0, 1.0);

            auto* leg = new TLegend(0.13, 0.79, 0.97, 0.925);
            leg->SetNColumns(2);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.032);
            for (int is = is_lo; is < is_hi; ++is) {
                ratios[is]->SetMarkerStyle(smark[is % smark.size()]);
                ratios[is]->SetMarkerColor(scol[is % scol.size()]);
                ratios[is]->SetLineColor(scol[is % scol.size()]);
                ratios[is]->SetLineWidth(2);
                leg->AddEntry(ratios[is], slabels[is].c_str(), "lp");
            }
            leg->Draw();
            std::vector<std::string> offscale;
            for (int is = is_lo; is < is_hi; ++is) {
                ratios[is]->Draw("E1 same");
                for (int i = 1; i <= ratios[is]->GetNbinsX(); ++i) {
                    const double v = ratios[is]->GetBinContent(i);
                    if (v <= ymax) continue;
                    const double x = ratios[is]->GetBinCenter(i);
                    auto* ar = new TArrow(x, ymax * 0.88, x, ymax * 0.985, 0.012, "|>");
                    ar->SetLineColor(scol[is % scol.size()]);
                    ar->SetFillColor(scol[is % scol.size()]);
                    ar->SetLineWidth(2);
                    ar->Draw();
                    offscale.push_back(Form("#DeltaR=%.2f: %.1f #pm %.1f", x, v,
                                            ratios[is]->GetBinError(i)));
                }
            }
            DrawHeadline(headline, 0.12, 0.965);
            if (!offscale.empty()) {
                TLatex note;
                note.SetNDC(); note.SetTextFont(42); note.SetTextSize(0.026);
                note.SetTextColor(kGray + 3);
                double y = 0.72;
                for (size_t i = 0; i < offscale.size(); i += 2) {
                    std::string txt = (i == 0) ? "above scale (arrows): " : "  ";
                    txt += offscale[i];
                    if (i + 1 < offscale.size()) txt += ", " + offscale[i + 1];
                    note.DrawLatex(0.45, y, txt.c_str());
                    y -= 0.035;
                }
            }
            SaveCanvas(c, dir3 + Form("step3_eps_dr_zoom_pair_pt_slices_%s.png",
                                      half == 0 ? "lowpt" : "highpt"));
        }
    }

    // ================================================================
    // Step 3 pair-eta dependence (round-5 #4): eps_dR(dR) in (pair pT, pair eta) cells.
    //   - two panel plots (zoom + full dR): each subplot = a pair-eta bin, each line = a
    //     pair-pT bin.
    //   - two tables: the large-dR plateau per (pair pT, pair eta) cell, and its fluctuation
    //     (weighted-mean stat error + weighted RMS scatter of the plateau-window bins).
    // Motivation: a plateau != 1 is manually normalized to 1 and the deviation is absorbed
    // into a systematic; quantify the deviation and its pair-pT / pair-eta dependence to
    // validate that and size the systematic.
    // ================================================================
    std::cout << "\n===== Step 3 pair-eta dependence (" << sample << ", " << wp_text << ") =====\n";
    {
        TH3D* h3zn = GetObj<TH3D>(fmc3, "h_mc_dr_zoom_vs_pt_eta_num");
        TH3D* h3zd = GetObj<TH3D>(fmc3, "h_mc_dr_zoom_vs_pt_eta_denom");
        TH3D* h3fn = GetObj<TH3D>(fmc3, "h_mc_dr_full_vs_pt_eta_num");
        TH3D* h3fd = GetObj<TH3D>(fmc3, "h_mc_dr_full_vs_pt_eta_denom");
        TH3D* h3zA = GetObj<TH3D>(fmc3, "h_mc_dr_zoom_vs_pt_eta_errA");
        TH3D* h3zB = GetObj<TH3D>(fmc3, "h_mc_dr_zoom_vs_pt_eta_errB");
        TH3D* h3fA = GetObj<TH3D>(fmc3, "h_mc_dr_full_vs_pt_eta_errA");
        TH3D* h3fB = GetObj<TH3D>(fmc3, "h_mc_dr_full_vs_pt_eta_errB");

        const int npt  = h3fn->GetYaxis()->GetNbins();  // coarse pair-pT (crossx)
        const int neta = h3fn->GetZaxis()->GetNbins();  // coarse pair-eta (crossx)

        auto pt_label  = [&](int iy){ return std::string(Form("%.1f < p_{T}^{pair} < %.1f GeV",
            h3fn->GetYaxis()->GetBinLowEdge(iy), h3fn->GetYaxis()->GetBinUpEdge(iy))); };
        auto eta_label = [&](int iz){ return std::string(Form("%.1f < #eta^{pair} < %.1f",
            h3fn->GetZaxis()->GetBinLowEdge(iz), h3fn->GetZaxis()->GetBinUpEdge(iz))); };

        // eps_dR(dR) for one (pt bin iy, eta bin iz) cell. Central value = num/denom;
        // ERROR = conditional/binomial-correct form (SetConditionalRatioErrors).
        // iz == 0 => integrate over ALL pair-eta (the pair-pT slices view), exactly as Step 4.
        auto cell_ratio = [](TH3D* hn, TH3D* hd, TH3D* ha, TH3D* hb,
                             int iy, int iz, int neta_all, const char* nm) -> TH1D* {
            const int zlo = (iz == 0) ? 1 : iz, zhi = (iz == 0) ? neta_all : iz;
            TH1D* n = hn->ProjectionX(Form("%s_n", nm), iy, iy, zlo, zhi, "e");
            TH1D* d = hd->ProjectionX(Form("%s_d", nm), iy, iy, zlo, zhi, "e");
            TH1D* a = ha->ProjectionX(Form("%s_a", nm), iy, iy, zlo, zhi, "e");
            TH1D* b = hb->ProjectionX(Form("%s_b", nm), iy, iy, zlo, zhi, "e");
            auto* r = (TH1D*)n->Clone(nm);
            r->SetDirectory(nullptr);
            r->Divide(d);
            SetConditionalRatioErrors(r, d, a, b);
            delete n; delete d; delete a; delete b;
            return r;
        };

        // 8 pair-pT bins: a 4-entry palette made bins 1&5, 2&6, 3&7, 4&8 indistinguishable.
        const std::vector<Color_t> ptcol  = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta,
                                             kOrange + 7, kCyan + 2, kViolet + 1, kBlack};
        const std::vector<Style_t> ptmark = {20, 21, 22, 23, 33, 34, 29, 24};

        // ---- two panel plots (zoom + full dR) ----
        struct Rng { TH3D* n; TH3D* d; TH3D* a; TH3D* b; std::string tag; double xhi; };
        const std::vector<Rng> rngs = {{h3zn, h3zd, h3zA, h3zB, "zoom", 1.0},
                                       {h3fn, h3fd, h3fA, h3fB, "full", 5.75}};
        // subplot grid: nrows >= ncols, nrows ~ sqrt(neta) (feedback_subplot_layout)
        const int ncol = (int)std::ceil(std::sqrt((double)neta));
        const int nrow = (int)std::ceil((double)neta / ncol);
        for (const auto& R : rngs) {
            // ONE y range for all nine panels of a canvas, built in a FIRST PASS over every
            // (pair pT, pair eta) cell -- exactly as the eps_dR distribution canvases below do.
            // Panel-by-panel autoscaling gave visible maxima of 2.4 / 2.0 / 1.8 / 1.6 on one
            // canvas, so two neighbouring panels of the same quantity were read on different
            // scales, and the two figure sets of that quantity followed opposite rules.
            // Range from CENTRAL values only, as those canvases do: the last wide-dR bins carry
            // errors of order 1 and including them stretches the axis to the 3.0 cap and flattens
            // the structure the figure exists for. Points above the cap keep their arrow AND are
            // listed on the canvas.
            std::vector<std::vector<TH1D*>> rs_all(neta + 1);
            double ymax = 0.;
            for (int iz = 1; iz <= neta; ++iz) {
                for (int iy = 1; iy <= npt; ++iy) {
                    TH1D* r = cell_ratio(R.n, R.d, R.a, R.b, iy, iz, neta,
                                         Form("s3pe_%s_%s_%d_%d", sample.c_str(), R.tag.c_str(), iy, iz));
                    rs_all[iz].push_back(r);
                    for (int i = 1; i <= r->GetNbinsX(); ++i)
                        ymax = std::max(ymax, r->GetBinContent(i));
                }
            }
            ymax = std::min(std::max(1.15, 1.15 * ymax), 3.0);   // cap; off-scale points arrowed

            // A RESERVED HEADER STRIP, with the grid of panels in a TPad below it -- the same
            // construction the eps_dR distribution canvases use. The eight pair-pT entries used
            // to be squeezed into ONE row across the full canvas width (SetNColumns(8) ~ 0.12 NDC
            // ~ 180 px per entry for a label like "104.0 < p_{T}^{pair} < 150.0 GeV"), so every
            // entry's text ran into the next entry's marker, the last one was clipped at the right
            // edge, and the only key identifying the eight curves was unreadable. Two rows of four
            // give each entry ~360 px.
            const int    kHeaderPx = 128;
            const int    canv_h    = 450 * nrow + kHeaderPx;
            const double hfrac     = (double)kHeaderPx / canv_h;
            TCanvas c(("c_s3_pteta_" + R.tag).c_str(), "", 500 * ncol, canv_h);
            auto* grid = new TPad(("grid_s3pe_" + R.tag).c_str(), "", 0., 0., 1., 1. - hfrac);
            grid->SetFillStyle(0); grid->Draw(); grid->cd(); grid->Divide(ncol, nrow);
            std::vector<TH1D*> rs_leg;   // series of the first panel, used for the canvas legend
            std::vector<std::string> offscale;
            for (int iz = 1; iz <= neta; ++iz) {
                grid->cd(iz);
                gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
                // Reserved strip ABOVE the frame for the pair-eta label, as on the eps_dR
                // distribution canvases: inside the frame the label sat at the same height as the
                // off-scale arrows, which are drawn just under the frame top, and the arrow was
                // painted straight through "-2.4 < eta^pair < -2.0".
                gPad->SetTopMargin(0.10);
                std::vector<TH1D*>& rs = rs_all[iz];
                DrawEffFrame(0.0, R.xhi, "#DeltaR", 0.0, ymax, cfg.eps_dr_text);
                DrawUnityLine(0.0, R.xhi);
                for (int iy = 0; iy < npt; ++iy) {
                    rs[iy]->SetMarkerStyle(ptmark[iy % ptmark.size()]);
                    rs[iy]->SetMarkerColor(ptcol[iy % ptcol.size()]);
                    rs[iy]->SetLineColor(ptcol[iy % ptcol.size()]);
                    rs[iy]->SetLineWidth(2);
                }
                if (iz == 1) rs_leg = rs;    // legend goes in the canvas top strip, see below
                for (int iy = 0; iy < npt; ++iy) {
                    TH1D* r = rs[iy];
                    r->Draw("E1 same");
                    for (int i = 1; i <= r->GetNbinsX(); ++i) {
                        if (r->GetBinContent(i) <= ymax) continue;
                        auto* ar = new TArrow(r->GetBinCenter(i), ymax * 0.88,
                                              r->GetBinCenter(i), ymax * 0.985, 0.008, "|>");
                        ar->SetLineColor(ptcol[iy % ptcol.size()]);
                        ar->SetFillColor(ptcol[iy % ptcol.size()]);
                        ar->Draw();
                        offscale.push_back(Form(
                            "#eta^{pair} #in [%.1f,%.1f), p_{T}^{pair} #in [%.1f,%.1f) GeV,"
                            " #DeltaR = %.2f: %.1f",
                            h3fn->GetZaxis()->GetBinLowEdge(iz), h3fn->GetZaxis()->GetBinUpEdge(iz),
                            h3fn->GetYaxis()->GetBinLowEdge(iy + 1),
                            h3fn->GetYaxis()->GetBinUpEdge(iy + 1),
                            r->GetBinCenter(i), r->GetBinContent(i)));
                    }
                }
                // eta-bin label in the reserved strip above the frame. It used to sit at y=0.94
                // of a pad with no top margin, i.e. under the canvas super-title; the strip keeps
                // it clear of BOTH the super-title (now in its own header) and the arrows inside
                // the frame.
                TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.050);
                tl.DrawLatex(0.15, 0.945, eta_label(iz).c_str());
            }
            c.cd(0);
            // Header strip, in PIXELS of this canvas so the layout does not depend on the grid
            // size. Row 1 = super-title, rows 2-3 = the two-row legend, row 4 = the off-scale
            // record. No `R.tag` ("zoom"/"full") on the canvas: which DeltaR range is shown is
            // visible on the axis, and the token is an internal name for two output files.
            auto ny = [&](double px) { return 1.0 - px / canv_h; };
            TLatex st; st.SetNDC(); st.SetTextFont(42); st.SetTextSize(22.0 / canv_h);
            st.DrawLatex(0.02, ny(26), (headline + ",  " + cfg.eps_dr_text +
                         " in (p_{T}^{pair}, #eta^{pair}) cells").c_str());
            if (!rs_leg.empty()) {
                // ONE legend for the whole canvas, in its own strip: in a panel the series cover
                // the frame, so an in-panel legend always sits on data.
                auto* leg = new TLegend(0.02, ny(100), 0.98, ny(38));
                leg->SetNColumns(4);
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(16.0 / canv_h);
                for (int iy = 0; iy < npt; ++iy)
                    leg->AddEntry(rs_leg[iy], pt_label(iy + 1).c_str(), "lp");
                leg->Draw();
            }
            if (!offscale.empty()) {
                // Terse record of every point the axis cap pushed off frame -- an arrow says a
                // point is missing, it does not say how far off it was.
                std::string note = "above the axis range: " + offscale[0];
                for (size_t k = 1; k < offscale.size() && k < 2; ++k) note += ";  " + offscale[k];
                if (offscale.size() > 2)
                    note += Form(";  ... (%d in total)", (int)offscale.size());
                TLatex os_; os_.SetNDC(); os_.SetTextFont(42); os_.SetTextSize(15.0 / canv_h);
                os_.SetTextColor(kGray + 3);
                os_.DrawLatex(0.02, ny(118), note.c_str());
            }
            SaveCanvas(c, dir3 + "step3_eps_dr_" + R.tag + "_pair_eta_pt.png");
        }

        // ---- tables (full-dR 3D, plateau window [kPlateauLo, kPlateauHi]) ----
        // Each cell is measured TWICE -- once in the nominal window and once in the retired
        // [1,4] one -- from the SAME ratio histogram, so the difference is purely the window.
        using Plat = PlateauCell;   // shared with the ROOT-file writer below
        auto cell_plateau = [&](int iy, int iz) -> Plat {
            TH1D* r = cell_ratio(h3fn, h3fd, h3fA, h3fB, iy, iz, neta,
                                 Form("plat_%s_%d_%d", sample.c_str(), iy, iz));
            const Plat p = PlateauFromRatio(r);
            delete r;
            return p;
        };

        std::vector<std::vector<Plat>> P(npt, std::vector<Plat>(neta));
        for (int iy = 1; iy <= npt; ++iy)
            for (int iz = 1; iz <= neta; ++iz) P[iy - 1][iz - 1] = cell_plateau(iy, iz);

        // ============================================================================
        // ROUND 9 (user request): the eps_dR DISTRIBUTIONS, one PNG per pair-pT bin, one
        // subplot per pair-eta bin, with THAT CELL'S plateau -- the very value the curve is
        // divided by before it is fitted -- drawn as a horizontal line. Three dR ranges, each
        // in its own subdirectory, because a single x range cannot show both the small-dR rise
        // (which lives below dR ~ 0.3) and the plateau region (dR in [2, 3.5]) legibly.
        //   0 to 1   : the fit domain, on the fine 0.05-wide dR bins
        //   0 to 2   : fine bins below 1, the wide bins above -- the same concatenation the fit
        //              plots use, so the two figure sets show literally the same points
        //   full     : out to the end of the wide-bin axis, so the plateau line sits ON the data
        // Superseded here: one canvas carrying all 8 pair-pT curves at once, which is what made
        // these unreadable.
        // ============================================================================
        {
            struct DrRange { std::string dir; double xhi; bool zoom_part; bool wide_part; };
            const std::vector<DrRange> dr_views = {
                {"dr_0_to_1",     1.0,  true,  false},
                {"dr_0_to_2",     2.0,  true,  true },
                {"dr_full_range", h3fn->GetXaxis()->GetXmax(), false, true }};

            // One graph per cell for a given view: fine bins below 1 and/or wide bins above.
            auto cell_graph = [&](const DrRange& R, int iy, int iz) -> TGraphErrors* {
                auto* g = new TGraphErrors();
                int k = 0;
                if (R.zoom_part) {
                    TH1D* rz = cell_ratio(h3zn, h3zd, h3zA, h3zB, iy, iz, neta,
                                          Form("d9z_%s_%d_%d", sample.c_str(), iy, iz));
                    for (int i = 1; i <= rz->GetNbinsX(); ++i) {
                        const double x = rz->GetBinCenter(i);
                        if (x > R.xhi || rz->GetBinContent(i) == 0.) continue;
                        g->SetPoint(k, x, rz->GetBinContent(i));
                        g->SetPointError(k, 0., rz->GetBinError(i));
                        ++k;
                    }
                    delete rz;
                }
                if (R.wide_part) {
                    TH1D* rf = cell_ratio(h3fn, h3fd, h3fA, h3fB, iy, iz, neta,
                                          Form("d9f_%s_%d_%d", sample.c_str(), iy, iz));
                    // below 1.0 the fine bins already cover the range -- never draw both
                    const double xlo_wide = R.zoom_part ? 1.0 : 0.0;
                    for (int i = 1; i <= rf->GetNbinsX(); ++i) {
                        const double x = rf->GetBinCenter(i);
                        if (x <= xlo_wide || x > R.xhi || rf->GetBinContent(i) == 0.) continue;
                        g->SetPoint(k, x, rf->GetBinContent(i));
                        g->SetPointError(k, 0., rf->GetBinError(i));
                        ++k;
                    }
                    delete rf;
                }
                return g;
            };

            const int ncol9 = (int)std::ceil(std::sqrt((double)neta));
            const int nrow9 = (int)std::ceil((double)neta / ncol9);
            int n_no_plateau = 0;
            for (const auto& R : dr_views) {
                const std::string vdir = dir3 + R.dir + "/";
                gSystem->mkdir(vdir.c_str(), kTRUE);
                for (int iy = 1; iy <= npt; ++iy) {
                    // y range shared by all panels of ONE png, from the drawn points and the
                    // plateau lines, so the nine pair-eta panels are directly comparable.
                    // Range from CENTRAL values only, as the fit canvases do: the last few
                    // wide-dR bins carry errors of order 1, and including them would stretch
                    // the axis to the 3.0 cap and flatten the structure this figure exists for.
                    double ylo = 1., yhi = 1.;
                    std::vector<TGraphErrors*> gs(neta + 1, nullptr);
                    for (int iz = 1; iz <= neta; ++iz) {
                        gs[iz] = cell_graph(R, iy, iz);
                        for (int i = 0; i < gs[iz]->GetN(); ++i) {
                            double x, y; gs[iz]->GetPoint(i, x, y);
                            ylo = std::min(ylo, y); yhi = std::max(yhi, y);
                        }
                        const Plat& pc = P[iy - 1][iz - 1];
                        if (pc.nb > 0) { ylo = std::min(ylo, pc.mean); yhi = std::max(yhi, pc.mean); }
                    }
                    const double span = std::max(yhi - ylo, 0.10);
                    ylo = std::max(0.0, ylo - 0.08 * span);
                    yhi = std::min(3.0, yhi + 0.22 * span);

                    // TWO header rows. The PbPb headline ("Pythia8 + HIJING overlay, Pb+Pb
                    // sqrt(s_NN) = 5.36 TeV, 0-5% (2023 conditions), Tight muons, ...") is far
                    // longer than the pp one, and a single-row header put it straight through the
                    // legend. Row 1 = headline, row 2 = the definition of the plotted quantity on
                    // the left and the legend on the right.
                    const int kHeaderPx = 118;
                    TCanvas c(Form("c_s3d9_%s_%d", R.dir.c_str(), iy), "",
                              520 * ncol9, 470 * nrow9 + kHeaderPx);
                    const double hfrac = (double)kHeaderPx / (470.0 * nrow9 + kHeaderPx);
                    auto* grid = new TPad(Form("grid_s3d9_%d", iy), "", 0., 0., 1., 1. - hfrac);
                    grid->SetFillStyle(0); grid->Draw(); grid->cd(); grid->Divide(ncol9, nrow9);

                    TLine* leg_line = nullptr;
                    std::vector<std::string> offscale;   // points pushed off the top of the frame
                    for (int iz = 1; iz <= neta; ++iz) {
                        grid->cd(iz);
                        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
                        gPad->SetTopMargin(0.10);
                        DrawEffFrame(0.0, R.xhi, "#DeltaR(#mu_{1}, #mu_{2})", ylo, yhi,
                                     cfg.eps_dr_text);
                        DrawUnityLine(0.0, R.xhi);

                        const Plat& pc = P[iy - 1][iz - 1];
                        TLine* pl = nullptr;
                        if (DrCorrPlateauUsable(pc.mean, pc.err)) {
                            // solid over the window where the window is on the canvas, dashed
                            // across the pad otherwise -- the line is a REFERENCE LEVEL there,
                            // not a measurement of the bins under it.
                            const bool in_range = (kPlateauLo < R.xhi);
                            pl = in_range ? new TLine(kPlateauLo, pc.mean,
                                                      std::min(kPlateauHi, R.xhi), pc.mean)
                                          : new TLine(0.0, pc.mean, R.xhi, pc.mean);
                            pl->SetLineColor(kBlue + 1);
                            pl->SetLineWidth(in_range ? 3 : 2);
                            if (!in_range) pl->SetLineStyle(2);
                            pl->Draw("same");
                            if (!leg_line) leg_line = pl;
                        } else if (R.dir == dr_views.front().dir) {   // count each cell ONCE, in the first view
                            ++n_no_plateau;
                        }

                        gs[iz]->SetMarkerStyle(20); gs[iz]->SetMarkerSize(0.8);
                        gs[iz]->SetMarkerColor(kBlack); gs[iz]->SetLineColor(kBlack);
                        gs[iz]->Draw("PZ same");

                        // A point above the y cap must be RECORDED, never silently dropped: on the
                        // overlay several small-dR points reach 8-11 in a panel that is otherwise
                        // drawn normally, and those are exactly the points this figure exists to
                        // show. Arrow in the panel + a terse canvas-level list, as the fit
                        // canvases already do.
                        for (int i = 0; i < gs[iz]->GetN(); ++i) {
                            double x, y; gs[iz]->GetPoint(i, x, y);
                            if (y <= yhi) continue;
                            auto* ar = new TArrow(x, ylo + 0.88 * (yhi - ylo),
                                                  x, ylo + 0.985 * (yhi - ylo), 0.008, "|>");
                            ar->SetLineColor(kBlack); ar->SetFillColor(kBlack); ar->Draw();
                            offscale.push_back(Form("#eta^{pair} #in [%.1f,%.1f), #DeltaR = %.2f: %.2f",
                                                    h3fn->GetZaxis()->GetBinLowEdge(iz),
                                                    h3fn->GetZaxis()->GetBinUpEdge(iz), x, y));
                        }

                        // Cell label AND plateau value both go ABOVE the frame, on one line. In
                        // the overlay the points and the off-scale arrows fill the frame edge to
                        // edge, so there is no free band inside it for an annotation -- and a
                        // number printed over the data is worse than one printed outside it.
                        TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.040);
                        tl.DrawLatex(0.15, 0.945,
                            (eta_label(iz) + ",   " +
                             (DrCorrPlateauUsable(pc.mean, pc.err)
                                  ? std::string(Form("plateau = %.4f #pm %.4f", pc.mean, pc.err))
                                  : std::string("plateau not measurable"))).c_str());
                    }
                    c.cd(0);
                    // Row 1: sample identity + the pair-pT cell. No eps symbol here -- every y
                    // axis on the canvas already carries it.
                    TLatex st; st.SetNDC(); st.SetTextFont(42); st.SetTextSize(0.016);
                    st.DrawLatex(0.03, 1. - 0.30 * hfrac,
                        (headline + ",  " + pt_label(iy)).c_str());
                    // Row 2 left: what the symbol on the y axis MEANS. The quantity is specific to
                    // this analysis, so the figure has to define it; §3.3 wording.
                    TLatex df; df.SetNDC(); df.SetTextFont(42); df.SetTextSize(0.015);
                    df.DrawLatex(0.03, 1. - 0.66 * hfrac,
                        (cfg.eps_dr_text + " = P(both #mu fire | #DeltaR) / "
                         "(#varepsilon_{1}#varepsilon_{2})").c_str());
                    if (!offscale.empty()) {
                        // Terse record of every point the axis cap pushed off frame.
                        std::string note = "above the axis range: " + offscale[0];
                        for (size_t k = 1; k < offscale.size() && k < 3; ++k) note += ";  " + offscale[k];
                        if (offscale.size() > 3)
                            note += Form(";  ... (%d in total)", (int)offscale.size());
                        TLatex os_; os_.SetNDC(); os_.SetTextFont(42); os_.SetTextSize(0.013);
                        os_.SetTextColor(kGray + 3);
                        os_.DrawLatex(0.03, 1. - 0.95 * hfrac, note.c_str());
                    }
                    if (!gs[1]->GetN() && !leg_line) { /* nothing drawn -- no legend to make */ }
                    else {
                        auto* lg = new TLegend(0.70, 1. - 0.90 * hfrac, 0.99, 1. - 0.42 * hfrac);
                        lg->SetNColumns(1);
                        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.015);
                        lg->AddEntry(gs[1], "measurement", "lp");
                        if (leg_line)
                            lg->AddEntry(leg_line,
                                Form("plateau, #DeltaR #in [%g, %g]", kPlateauLo, kPlateauHi), "l");
                        lg->Draw();
                    }
                    // ONE DECIMAL: at %.0f the log edges 11.54 / 16.65 / 24.01 print as 12/17/24
                    // while the canvas label says 11.5/16.6/24.0 -- a name that contradicts the
                    // figure is the binning-drift the repo rule exists to prevent.
                    SaveCanvas(c, vdir + Form("step3_eps_dr_pairpt_%.1f_%.1f.png",
                                              h3fn->GetYaxis()->GetBinLowEdge(iy),
                                              h3fn->GetYaxis()->GetBinUpEdge(iy)));
                    for (int iz = 1; iz <= neta; ++iz) delete gs[iz];
                }
            }
            printf("  Step-3 dR-correction distributions: %d x %d PNGs written under %s"
                   " (dr_0_to_1 / dr_0_to_2 / dr_full_range); %d of %d cells have no measurable"
                   " plateau and carry no plateau line\n",
                   (int)dr_views.size(), npt, dir3.c_str(), n_no_plateau, npt * neta);
        }

        // Table A: plateau value +- stat error ; Table B: fluctuation (stat err + RMS scatter)
        auto eta_hdr = [&](int iz){ return Form("[%.1f,%.1f)",
            h3fn->GetZaxis()->GetBinLowEdge(iz), h3fn->GetZaxis()->GetBinUpEdge(iz)); };
        auto pt_hdr  = [&](int iy){ return Form("pTpair[%.1f,%.1f)",
            h3fn->GetYaxis()->GetBinLowEdge(iy), h3fn->GetYaxis()->GetBinUpEdge(iy)); };

        auto write_table_A = [&](std::ostream& os){
            os << "# Step-3 large-dR plateau eps_dR (weighted mean over dR in ["
               << kPlateauLo << "," << kPlateauHi << "]) per (pair pT, pair eta) cell.\n";
            os << "# sample=" << sample << "  WP=" << wp_text << "  " << cfg.sample_text << "\n";
            os << "# value = plateau +- stat_error   (a value != 1 is the manual-normalization"
                  " correction; the plateau is normalized to 1 and |value-1| feeds the systematic)\n";
            os << std::left << std::setw(20) << "pair-eta \\ pair-pT";
            for (int iy = 1; iy <= npt; ++iy) os << std::setw(20) << pt_hdr(iy);
            os << "\n";
            for (int iz = 1; iz <= neta; ++iz) {
                os << std::left << std::setw(20) << eta_hdr(iz);
                for (int iy = 1; iy <= npt; ++iy) {
                    const Plat& p = P[iy - 1][iz - 1];
                    os << std::setw(20) << (p.nb > 0 ? Form("%.4f+-%.4f", p.mean, p.err) : "--");
                }
                os << "\n";
            }
        };
        auto write_table_B = [&](std::ostream& os){
            os << "# Step-3 plateau FLUCTUATION per (pair pT, pair eta) cell (dR in ["
               << kPlateauLo << "," << kPlateauHi << "]).\n";
            os << "# sample=" << sample << "  WP=" << wp_text << "\n";
            os << "# each cell = stat_err(mean) / rms_scatter (weighted RMS of the plateau-window"
                  " bins about their mean) / n_bins / window_syst.\n";
            os << "# rms_scatter measures how flat/noisy the plateau is; compare |mean-1| to it +"
                  " stat_err to judge whether the deviation is significant (systematic size).\n";
            os << "# window_syst = |plateau[" << kPlateauLo << "," << kPlateauHi << "] - plateau["
               << kPlateauSystLo << "," << kPlateauSystHi << "]|, the plateau-window"
                  " normalization systematic (-- = not evaluable).\n";
            os << std::left << std::setw(28) << "pair-eta \\ pair-pT";
            for (int iy = 1; iy <= npt; ++iy) os << std::setw(28) << pt_hdr(iy);
            os << "\n";
            for (int iz = 1; iz <= neta; ++iz) {
                os << std::left << std::setw(28) << eta_hdr(iz);
                for (int iy = 1; iy <= npt; ++iy) {
                    const Plat& p = P[iy - 1][iz - 1];
                    os << std::setw(28) << (p.nb > 0
                        ? Form("%.4f/%.4f/%d/%s", p.err, p.rms, p.nb,
                               p.syst >= 0. ? Form("%.4f", p.syst) : "--")
                        : "--");
                }
                os << "\n";
            }
        };
        write_table_A(std::cout);
        write_table_B(std::cout);
        { std::ofstream ofa(dir3 + "step3_plateau_pair_eta_pt.txt");           write_table_A(ofa); }
        { std::ofstream ofb(dir3 + "step3_plateau_fluctuation_pair_eta_pt.txt"); write_table_B(ofb); }
        std::cout << "  wrote " << dir3 << "step3_plateau_pair_eta_pt.txt + "
                  << "step3_plateau_fluctuation_pair_eta_pt.txt\n";

        // Same numbers, machine-readable, for fit_dr_corrections.cxx. RECREATE: Step 3 always
        // runs, so it owns the file and Step 4 appends to it below.
        WritePlateauRootFile(DrCorrPlateauFile(id, use_tight_wp), /*recreate=*/true, 3, h3fn, P,
                             plateau.first, plateau.second, plateau_syst,
                             sample, wp_text, "eps_dR");

        // ---- ROUND 9: the SAME-SIGN and OPPOSITE-SIGN series -----------------------------
        // Each sign gets its OWN plateau, measured exactly as the sign-integrated one and from
        // its own histograms -- normalizing a same-sign curve by an opposite-sign-dominated
        // plateau would import the very charge dependence the split exists to test. Written
        // into the same plateau file under sign-tagged keys; a sample filled before round 9
        // simply has no per-sign histograms and is skipped with a note.
        for (const std::string& sgn : {std::string("ss"), std::string("os")}) {
            const std::string hp = "h_mc_dr_" + sgn + "_";
            TH3D* sn = GetObjOrNull<TH3D>(fmc3, hp + "full_vs_pt_eta_num");
            TH3D* sd = GetObjOrNull<TH3D>(fmc3, hp + "full_vs_pt_eta_denom");
            TH3D* sA = GetObjOrNull<TH3D>(fmc3, hp + "full_vs_pt_eta_errA");
            TH3D* sB = GetObjOrNull<TH3D>(fmc3, hp + "full_vs_pt_eta_errB");
            TH1D* i_n = GetObjOrNull<TH1D>(fmc3, hp + "full_num");
            TH1D* i_d = GetObjOrNull<TH1D>(fmc3, hp + "full_denom");
            TH1D* i_A = GetObjOrNull<TH1D>(fmc3, hp + "full_errA");
            TH1D* i_B = GetObjOrNull<TH1D>(fmc3, hp + "full_errB");
            if (!sn || !sd || !sA || !sB || !i_n || !i_d || !i_A || !i_B) {
                std::cout << "  [Step 3] no " << hp << "* histograms -- sign-separated plateaus "
                          << "skipped for this sample (refill to produce them)\n";
                continue;
            }
            auto* ri = (TH1D*)i_n->Clone(("s3incl_" + sgn).c_str());
            ri->SetDirectory(nullptr);
            ri->Divide(i_d);
            SetConditionalRatioErrors(ri, i_d, i_A, i_B);
            const auto pi  = PlateauWeightedMean(ri, kPlateauLo, kPlateauHi);
            const auto pia = PlateauWeightedMean(ri, kPlateauSystLo, kPlateauSystHi);
            delete ri;

            std::vector<std::vector<Plat>> Ps(npt, std::vector<Plat>(neta));
            for (int iy = 1; iy <= npt; ++iy)
                for (int iz = 1; iz <= neta; ++iz) {
                    TH1D* r = cell_ratio(sn, sd, sA, sB, iy, iz, neta,
                                         Form("plat3%s_%s_%d_%d", sgn.c_str(), sample.c_str(), iy, iz));
                    Ps[iy - 1][iz - 1] = PlateauFromRatio(r);
                    delete r;
                }
            WritePlateauRootFile(DrCorrPlateauFile(id, use_tight_wp), /*recreate=*/false, 3, sn, Ps,
                                 pi.first, pi.second, std::fabs(pi.first - pia.first),
                                 sample, wp_text, "eps_dR", sgn);
            printf("  %s (%s pairs) inclusive plateau = %.4f +- %.4f\n", sample.c_str(),
                   sgn == "ss" ? "same-sign" : "opposite-sign", pi.first, pi.second);
        }
    }

    // ================================================================
    // RAW joint trigger probability per pair charge combination (round 9)
    // ================================================================
    // WHY THIS FIGURE IS PART OF THE COMMITTED CHAIN. The Step-3 dR correction splits by pair
    // charge below dR ~ 0.15 (mc_trigger_efficiency.md R25). Two explanations competed: a genuine
    // L1 close-by-RoI effect, or an eps_MC mis-parameterisation that SQUARES in the 1/(eps1 eps2)
    // numerator weight for a close SAME-sign pair (q.eta_1 ~ q.eta_2) while partly CANCELLING for
    // a close OPPOSITE-sign pair (q.eta_1 ~ -q.eta_2). The discriminator is the numerator that
    // carries NO eps_MC at all -- the SAME trigger-passing node weighted by the plain MC weight
    // (booked as "numraw" in FillMCTrigEffHists.cxx) -- because the second explanation lives
    // entirely in the weight. R25's refutation of it rests on this quantity and on nothing else,
    // so the chain must DRAW it; until it did, the result was not reproducible by running the
    // committed code. Graceful skip for a sample filled before that numerator was booked.
    {
        TH1D* n_ss = GetObjOrNull<TH1D>(fmc3, "h_mc_dr_ss_zoom_numraw");
        TH1D* d_ss = GetObjOrNull<TH1D>(fmc3, "h_mc_dr_ss_zoom_denom");
        TH1D* n_os = GetObjOrNull<TH1D>(fmc3, "h_mc_dr_os_zoom_numraw");
        TH1D* d_os = GetObjOrNull<TH1D>(fmc3, "h_mc_dr_os_zoom_denom");
        if (!n_ss || !d_ss || !n_os || !d_os) {
            std::cout << "  [Step 3] the per-sign RAW joint trigger probability is NOT available "
                         "in this file (h_mc_dr_{ss,os}_zoom_numraw missing) -- figure skipped; "
                         "refill this sample to produce it\n";
        } else {
            // P = sum_fired w / sum_all w. With a_i = w_i (no inverse weight) the round-7
            // conditional variance A - R*B specialises to A = B = sum_fired w^2, which Sumw2
            // already stores as the numerator's bin error squared. Feeding that to the SHARED
            // SetConditionalRatioErrors keeps one implementation of the error formula instead of
            // a second, drifting copy.
            auto raw_prob = [](TH1D* n, TH1D* d, const char* nm) -> TH1D* {
                auto* ab = (TH1D*)n->Clone(Form("%s_ab", nm));
                ab->SetDirectory(nullptr);
                for (int i = 1; i <= ab->GetNbinsX(); ++i) {
                    const double e = n->GetBinError(i);
                    ab->SetBinContent(i, e * e);
                }
                auto* r = (TH1D*)n->Clone(nm);
                r->SetDirectory(nullptr);
                r->Divide(d);
                SetConditionalRatioErrors(r, d, ab, ab);
                delete ab;
                return r;
            };
            TH1D* p_ss = raw_prob(n_ss, d_ss, "praw_ss");
            TH1D* p_os = raw_prob(n_os, d_os, "praw_os");

            // Same-sign and opposite-sign pairs are DISJOINT samples, so the two probabilities
            // are statistically independent and their relative errors add in quadrature.
            auto* p_rat = (TH1D*)p_ss->Clone("praw_ss_over_os");
            p_rat->SetDirectory(nullptr);
            for (int i = 1; i <= p_rat->GetNbinsX(); ++i) {
                const double a = p_ss->GetBinContent(i), b = p_os->GetBinContent(i);
                const double ea = p_ss->GetBinError(i),  eb = p_os->GetBinError(i);
                if (!(a > 0.) || !(b > 0.)) { p_rat->SetBinContent(i, 0.); p_rat->SetBinError(i, 0.); continue; }
                const double v = a / b;
                p_rat->SetBinContent(i, v);
                p_rat->SetBinError(i, v * std::sqrt((ea / a) * (ea / a) + (eb / b) * (eb / b)));
            }

            // Points, not bin contents: an empty denominator bin is genuinely undefined and must
            // not be drawn at zero.
            auto to_graph = [](TH1D* h, TH1D* d) {
                auto* g = new TGraphErrors();
                int k = 0;
                for (int i = 1; i <= h->GetNbinsX(); ++i) {
                    if (d->GetBinContent(i) <= 0. || h->GetBinContent(i) <= 0.) continue;
                    g->SetPoint(k, h->GetBinCenter(i), h->GetBinContent(i));
                    g->SetPointError(k, 0., h->GetBinError(i));
                    ++k;
                }
                return g;
            };
            auto* g_ss  = to_graph(p_ss,  d_ss);
            auto* g_os  = to_graph(p_os,  d_os);
            auto* g_rat = to_graph(p_rat, d_ss);

            const double xlo = p_ss->GetXaxis()->GetXmin();
            const double xhi = p_ss->GetXaxis()->GetXmax();
            auto span_of = [](TGraphErrors* g, double& lo, double& hi) {
                for (int i = 0; i < g->GetN(); ++i) {
                    double x, y; g->GetPoint(i, x, y);
                    lo = std::min(lo, y); hi = std::max(hi, y);
                }
            };
            double ymin = 1e30, ymax = -1e30;
            span_of(g_ss, ymin, ymax); span_of(g_os, ymin, ymax);
            if (ymin > ymax) { ymin = 0.; ymax = 1.; }
            double rmin = 1e30, rmax = -1e30;
            span_of(g_rat, rmin, rmax);
            if (rmin > rmax) { rmin = 0.; rmax = 2.; }
            const double ys = std::max(ymax - ymin, 0.05);
            const double rs = std::max(rmax - rmin, 0.10);

            // The trigger requirement itself differs by beam: pp uses the dimuon 2mu4 decision,
            // the Pb+Pb overlay requires each muon to pass mu4 (FillMCTrigEffHists.cxx
            // `trig_cond`). Name the requirement the reader is actually looking at.
            const bool cross = (sample == "overlay");
            const std::string ytitle = cross ? "P(both #mu pass mu4 | #mu pair)"
                                             : "P(2mu4 | #mu pair)";

            TCanvas c("c_step3_raw_joint_prob", "", 900, 780);
            TPad body("p_s3raw_body", "", 0., 0., 1., 0.93);   // headline gets its own strip
            body.Draw();
            body.cd();
            auto pads = SplitPadForRatio("s3raw", false);

            pads.first->cd();
            // Reserved strip above the frame: the two series fill the pad, so an in-frame legend
            // would sit on the data (and a white backing would only hide points).
            gPad->SetTopMargin(0.20);
            DrawEffFrame(xlo, xhi, "", std::max(0., ymin - 0.15 * ys),
                         std::min(1.05, ymax + 0.15 * ys), ytitle, true);
            g_os->SetMarkerColor(kBlue + 1); g_os->SetLineColor(kBlue + 1);
            g_os->SetMarkerStyle(20); g_os->SetMarkerSize(1.0); g_os->SetLineWidth(2);
            g_ss->SetMarkerColor(kRed + 1);  g_ss->SetLineColor(kRed + 1);
            g_ss->SetMarkerStyle(21); g_ss->SetMarkerSize(1.0); g_ss->SetLineWidth(2);
            g_os->Draw("PZ same");
            g_ss->Draw("PZ same");

            TLatex eq;
            eq.SetNDC(); eq.SetTextFont(42); eq.SetTextSize(0.046);
            eq.DrawLatex(0.15, 0.925,
                         "P = #Sigma w_{MC} (pairs that fired) / #Sigma w_{MC} (all pairs)");
            auto* leg = new TLegend(0.15, 0.815, 0.95, 0.885);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetNColumns(2);
            leg->SetTextSize(0.044);
            leg->AddEntry(g_ss, "same sign, P_{SS}", "lp");
            leg->AddEntry(g_os, "opposite sign, P_{OS}", "lp");
            leg->Draw();

            pads.second->cd();
            DrawRatioFrame(xlo, xhi, "#DeltaR", "P_{SS} / P_{OS}",
                           std::max(0., rmin - 0.15 * rs), rmax + 0.15 * rs);
            g_rat->SetMarkerColor(kBlack); g_rat->SetLineColor(kBlack);
            g_rat->SetMarkerStyle(20); g_rat->SetMarkerSize(1.0); g_rat->SetLineWidth(2);
            g_rat->Draw("PZ same");

            c.cd();
            DrawHeadline(headline, 0.12, 0.962, 0.040);
            SaveCanvas(c, dir3 + "step3_raw_joint_trigger_probability_by_sign.png");

            // The numbers behind the figure, so the charge split can be checked against values.
            std::cout << "  RAW joint trigger probability (no eps_MC weighting), " << sample
                      << ", " << wp_text << ":\n";
            for (int i = 1; i <= p_ss->GetNbinsX(); ++i) {
                if (p_ss->GetBinContent(i) <= 0. || p_os->GetBinContent(i) <= 0.) continue;
                if (p_ss->GetBinCenter(i) > 0.5) break;
                printf("    dR = %.3f : same sign %.4f, opposite sign %.4f, ratio %.3f\n",
                       p_ss->GetBinCenter(i), p_ss->GetBinContent(i), p_os->GetBinContent(i),
                       p_rat->GetBinContent(i));
            }
            delete p_ss; delete p_os; delete p_rat;
        }
    }

    // ================================================================
    // Step 4 (§3.4): single-leg ΔR correction ε_single(ΔR) = inverse-weighted num / denom.
    //   Leg-level analog of Step 3. Deliverable for the PbPb UNION linear terms (overlay);
    //   for pp it is a VALIDATION only (2mu4 product absorbs the single-leg ΔR into ε_ΔR^2mu4)
    //   -- annotated on the pp panels. Graceful skip if _step4.root is absent.
    // ================================================================
    {
        const std::string step4_path =
            cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf
                                 + MCTrigEffPairPt::FileSuffix() + "_step4.root";
        TFile* fmc4 = TFile::Open(step4_path.c_str(), "READ");
        if (!fmc4 || fmc4->IsZombie()) {
            std::cout << "\n[Step 4] no " << step4_path << " -- skipping single-leg ΔR plots "
                      << "(run FillMCTrigEffHists do_step4=true first).\n";
        } else {
            std::cout << "\n===== Step 4 (" << sample << ", " << wp_text << ") =====\n";
            const std::string dir4 = cfg.out_base + "step4_dr_correction_singles/" + wp_dir;
            gSystem->mkdir(dir4.c_str(), kTRUE);
            const std::string eps_single_text = "#varepsilon_{#DeltaR}^{single}";
            // Used ONLY in the .txt tables below (a working document, where the downstream role
            // is exactly the context a reader needs) -- never on a canvas.
            const bool is_deliverable = (sample == "overlay");   // PbPb union; else validation
            // NO downstream-usage note on the canvas (user decision 2026-08-04): whether this
            // curve is applied to the PbPb union or kept as a pp validation is an analysis
            // decision recorded in mc_trigger_efficiency.md, not a property of the measurement
            // the figure shows.

            // Central value = num/denom; ERROR = conditional/binomial-correct form INCLUDING the
            // leg-leg covariance (both legs of a pair share a dR bin and their trigger decisions
            // are correlated -- that correlation is what Step 3 measures).
            auto MakeRatio4 = [&](const std::string& tag) -> TH1D* {
                TH1D* num = GetObj<TH1D>(fmc4, "h_mc_single_dr_" + tag + "_num");
                TH1D* den = GetObj<TH1D>(fmc4, "h_mc_single_dr_" + tag + "_denom");
                TH1D* eA  = GetObj<TH1D>(fmc4, "h_mc_single_dr_" + tag + "_errA");
                TH1D* eB  = GetObj<TH1D>(fmc4, "h_mc_single_dr_" + tag + "_errB");
                TH1D* cP  = GetObj<TH1D>(fmc4, "h_mc_single_dr_" + tag + "_covP");
                TH1D* cQ  = GetObj<TH1D>(fmc4, "h_mc_single_dr_" + tag + "_covQ");
                auto* r = (TH1D*)num->Clone(("r_single_dr_" + tag).c_str());
                r->SetDirectory(nullptr);
                r->Divide(den);
                SetConditionalRatioErrors(r, den, eA, eB, cP, cQ);
                return r;
            };
            TH1D* r4_zoom = MakeRatio4("zoom");
            TH1D* r4_full = MakeRatio4("full");

            const auto plat4 = PlateauWeightedMean(r4_full, kPlateauLo, kPlateauHi);
            const auto plat4_alt = PlateauWeightedMean(r4_full, kPlateauSystLo, kPlateauSystHi);
            const double plat4_syst = std::fabs(plat4.first - plat4_alt.first);
            printf("  %s single-leg plateau (weighted mean, dR in [%g,%g]): %.4f +- %.4f"
                   "   [window syst vs dR in [%g,%g] = %.4f]\n",
                   sample.c_str(), kPlateauLo, kPlateauHi, plat4.first, plat4.second,
                   kPlateauSystLo, kPlateauSystHi, plat4_syst);

            auto DrawStep4 = [&](TH1D* r, double xlo, double xhi, const std::string& png,
                                 bool plateau_in_range) {
                TCanvas c(("c_" + png).c_str(), "", 900, 700);
                gPad->SetLeftMargin(0.12);
                gPad->SetBottomMargin(0.12);
                double ymax = 0.;
                for (int i = 1; i <= r->GetNbinsX(); ++i)
                    ymax = std::max(ymax, r->GetBinContent(i) + r->GetBinError(i));
                ymax = std::max(1.15, 1.15 * ymax);
                DrawEffFrame(xlo, xhi, "#DeltaR", 0.0, ymax, eps_single_text);
                DrawUnityLine(xlo, xhi);
                TLine* plateau_line = nullptr;
                if (plateau_in_range) {
                    plateau_line = new TLine(kPlateauLo, plat4.first, kPlateauHi, plat4.first);
                    plateau_line->SetLineColor(kBlue + 1); plateau_line->SetLineWidth(3);
                } else {
                    plateau_line = new TLine(xlo, plat4.first, xhi, plat4.first);
                    plateau_line->SetLineColor(kBlue + 1); plateau_line->SetLineWidth(2);
                    plateau_line->SetLineStyle(3);
                }
                plateau_line->Draw("same");
                r->SetMarkerStyle(20);
                r->SetMarkerColor(kMCColor);
                r->SetLineColor(kMCColor);
                r->SetLineWidth(2);
                r->Draw("E1 same");
                auto* leg4 = new TLegend(0.42, 0.40, 0.90, 0.52);
                leg4->SetBorderSize(0); leg4->SetFillStyle(0); leg4->SetTextSize(0.033);
                leg4->AddEntry(r, "measurement", "lp");
                leg4->AddEntry(plateau_line, "plateau", "l");
                leg4->Draw();
                DrawHeadline(headline);
                // Same placement rule as Step 3: annotation in the empty lower-right quadrant.
                TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.035);
                tl.DrawLatex(0.42, 0.32, (eps_single_text +
                    " = P(#mu fires mu4 | #DeltaR) / #varepsilon(p_{T}, q#upoint#eta)").c_str());
                tl.DrawLatex(0.42, 0.25, Form("plateau #LT#DeltaR#in[%g,%g]#GT = %.4f #pm %.4f",
                                              kPlateauLo, kPlateauHi, plat4.first, plat4.second));
                SaveCanvas(c, dir4 + png + ".png");
            };
            DrawStep4(r4_zoom, 0.0, 1.0,  "step4_eps_dr_single_zoom", false);
            DrawStep4(r4_full, 0.0, 5.75, "step4_eps_dr_single_full", true);

            // ---- pair-eta x pair-pT breakdown from the 3D (plateau-stability systematic) ----
            TH3D* h3zn = GetObj<TH3D>(fmc4, "h_mc_single_dr_zoom_vs_pt_eta_num");
            TH3D* h3zd = GetObj<TH3D>(fmc4, "h_mc_single_dr_zoom_vs_pt_eta_denom");
            TH3D* h3fn = GetObj<TH3D>(fmc4, "h_mc_single_dr_full_vs_pt_eta_num");
            TH3D* h3fd = GetObj<TH3D>(fmc4, "h_mc_single_dr_full_vs_pt_eta_denom");
            TH3D* h3zA = GetObj<TH3D>(fmc4, "h_mc_single_dr_zoom_vs_pt_eta_errA");
            TH3D* h3zB = GetObj<TH3D>(fmc4, "h_mc_single_dr_zoom_vs_pt_eta_errB");
            TH3D* h3zP = GetObj<TH3D>(fmc4, "h_mc_single_dr_zoom_vs_pt_eta_covP");
            TH3D* h3zQ = GetObj<TH3D>(fmc4, "h_mc_single_dr_zoom_vs_pt_eta_covQ");
            TH3D* h3fA = GetObj<TH3D>(fmc4, "h_mc_single_dr_full_vs_pt_eta_errA");
            TH3D* h3fB = GetObj<TH3D>(fmc4, "h_mc_single_dr_full_vs_pt_eta_errB");
            TH3D* h3fP = GetObj<TH3D>(fmc4, "h_mc_single_dr_full_vs_pt_eta_covP");
            TH3D* h3fQ = GetObj<TH3D>(fmc4, "h_mc_single_dr_full_vs_pt_eta_covQ");
            const int npt  = h3fn->GetYaxis()->GetNbins();
            const int neta = h3fn->GetZaxis()->GetNbins();
            auto pt_label  = [&](int iy){ return std::string(Form("%.1f < p_{T}^{pair} < %.1f GeV",
                h3fn->GetYaxis()->GetBinLowEdge(iy), h3fn->GetYaxis()->GetBinUpEdge(iy))); };
            auto eta_label = [&](int iz){ return std::string(Form("%.1f < #eta^{pair} < %.1f",
                h3fn->GetZaxis()->GetBinLowEdge(iz), h3fn->GetZaxis()->GetBinUpEdge(iz))); };
            // eps_single(dR) for one cell; iz=0 => integrate over ALL eta (pair-pT slice)
            auto cell_ratio = [](TH3D* hn, TH3D* hd, TH3D* ha, TH3D* hb, TH3D* hp, TH3D* hq,
                                 int iy, int iz, int neta_all, const char* nm) -> TH1D* {
                const int zlo = (iz == 0) ? 1 : iz, zhi = (iz == 0) ? neta_all : iz;
                TH1D* n = hn->ProjectionX(Form("%s_n", nm), iy, iy, zlo, zhi, "e");
                TH1D* d = hd->ProjectionX(Form("%s_d", nm), iy, iy, zlo, zhi, "e");
                TH1D* a = ha->ProjectionX(Form("%s_a", nm), iy, iy, zlo, zhi, "e");
                TH1D* b = hb->ProjectionX(Form("%s_b", nm), iy, iy, zlo, zhi, "e");
                TH1D* p = hp->ProjectionX(Form("%s_p", nm), iy, iy, zlo, zhi, "e");
                TH1D* q = hq->ProjectionX(Form("%s_q", nm), iy, iy, zlo, zhi, "e");
                auto* r = (TH1D*)n->Clone(nm);
                r->SetDirectory(nullptr);
                r->Divide(d);
                SetConditionalRatioErrors(r, d, a, b, p, q);
                delete n; delete d; delete a; delete b; delete p; delete q;
                return r;
            };
            // 8 pair-pT bins since round 8 -- one distinct colour+marker per bin.
            const std::vector<Color_t> ptcol  = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta,
                                                 kOrange + 7, kCyan + 2, kViolet + 1, kBlack};
            const std::vector<Style_t> ptmark = {20, 21, 22, 23, 33, 34, 29, 24};

            // pair-pT slices (integrate eta): zoom dR
            {
                // Build every series ONCE, then draw the lower and upper halves of the pair-pT
                // bins into SEPARATE PNGs (round 8, user): 8 overlaid series in one pad are
                // unreadable. y-range computed per half so each file is scaled to its own data.
                std::vector<TH1D*> rs;
                for (int iy = 1; iy <= npt; ++iy)
                    rs.push_back(cell_ratio(h3zn, h3zd, h3zA, h3zB, h3zP, h3zQ, iy, 0, neta,
                                            Form("s4pt_%s_%d", sample.c_str(), iy)));
                const int nhalf4 = (npt + 1) / 2;
                for (int half = 0; half < 2; ++half) {
                    const int lo = half * nhalf4, hi = std::min(npt, (half + 1) * nhalf4);
                    if (lo >= hi) continue;
                    TCanvas c(Form("c_step4_ptslices_%d", half), "", 900, 700);
                    gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
                    gPad->SetTopMargin(0.22);   // reserved legend strip, as in Step 3
                    double ymax = 0.;
                    for (int iy = lo; iy < hi; ++iy)
                        for (int i = 1; i <= rs[iy]->GetNbinsX(); ++i)
                            ymax = std::max(ymax, rs[iy]->GetBinContent(i) + rs[iy]->GetBinError(i));
                    ymax = std::min(std::max(1.15, 1.15 * ymax), 3.0);
                    DrawEffFrame(0.0, 1.0, "#DeltaR", 0.0, ymax, eps_single_text);
                    DrawUnityLine(0.0, 1.0);
                    auto* leg = new TLegend(0.13, 0.79, 0.97, 0.925);
                    leg->SetNColumns(2);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->SetTextSize(0.032);
                    for (int iy = lo; iy < hi; ++iy) {
                        rs[iy]->SetMarkerStyle(ptmark[iy % ptmark.size()]);
                        rs[iy]->SetMarkerColor(ptcol[iy % ptcol.size()]);
                        rs[iy]->SetLineColor(ptcol[iy % ptcol.size()]);
                        rs[iy]->SetLineWidth(2);
                        leg->AddEntry(rs[iy], pt_label(iy + 1).c_str(), "lp");
                    }
                    leg->Draw();
                    std::vector<std::string> offscale;
                    for (int iy = lo; iy < hi; ++iy) {
                        rs[iy]->Draw("E1 same");
                        for (int i = 1; i <= rs[iy]->GetNbinsX(); ++i) {
                            if (rs[iy]->GetBinContent(i) <= ymax) continue;
                            auto* ar = new TArrow(rs[iy]->GetBinCenter(i), ymax * 0.88,
                                                  rs[iy]->GetBinCenter(i), ymax * 0.985, 0.012, "|>");
                            ar->SetLineColor(ptcol[iy % ptcol.size()]);
                            ar->SetFillColor(ptcol[iy % ptcol.size()]); ar->Draw();
                            offscale.push_back(Form("#DeltaR=%.2f: %.1f #pm %.1f", rs[iy]->GetBinCenter(i),
                                                    rs[iy]->GetBinContent(i), rs[iy]->GetBinError(i)));
                        }
                    }
                    DrawHeadline(headline, 0.12, 0.965);
                    if (!offscale.empty()) {
                        TLatex note; note.SetNDC(); note.SetTextFont(42); note.SetTextSize(0.026);
                        note.SetTextColor(kGray + 3);
                        double y = 0.72;   // inside the frame, below the legend strip
                        for (size_t i = 0; i < offscale.size(); i += 2) {
                            std::string txt = (i == 0) ? "above scale (arrows): " : "  ";
                            txt += offscale[i];
                            if (i + 1 < offscale.size()) txt += ", " + offscale[i + 1];
                            note.DrawLatex(0.45, y, txt.c_str()); y -= 0.035;
                        }
                    }
                    SaveCanvas(c, dir4 + Form("step4_eps_dr_single_zoom_pair_pt_slices_%s.png",
                                              half == 0 ? "lowpt" : "highpt"));
                }
            }

            // pair-eta panels (each subplot = eta bin, each line = pair-pT bin), zoom + full
            struct Rng { TH3D* n; TH3D* d; TH3D* a; TH3D* b; TH3D* p; TH3D* q; std::string tag; double xhi; };
            const std::vector<Rng> rngs = {{h3zn, h3zd, h3zA, h3zB, h3zP, h3zQ, "zoom", 1.0},
                                          {h3fn, h3fd, h3fA, h3fB, h3fP, h3fQ, "full", 5.75}};
            const int ncol = (int)std::ceil(std::sqrt((double)neta));
            const int nrow = (int)std::ceil((double)neta / ncol);
            for (const auto& R : rngs) {
                // ONE y range for all nine panels, from a FIRST PASS over every (pair pT, pair
                // eta) cell, and CENTRAL values only -- identical rule to the Step-3 panels and
                // to the eps_dR distribution canvases. Per-panel autoscaling made neighbouring
                // panels of the same quantity read on different scales.
                std::vector<std::vector<TH1D*>> rs_all(neta + 1);
                double ymax = 0.;
                for (int iz = 1; iz <= neta; ++iz) {
                    for (int iy = 1; iy <= npt; ++iy) {
                        TH1D* r = cell_ratio(R.n, R.d, R.a, R.b, R.p, R.q, iy, iz, neta,
                                             Form("s4pe_%s_%s_%d_%d", sample.c_str(), R.tag.c_str(), iy, iz));
                        rs_all[iz].push_back(r);
                        for (int i = 1; i <= r->GetNbinsX(); ++i)
                            ymax = std::max(ymax, r->GetBinContent(i));
                    }
                }
                ymax = std::min(std::max(1.15, 1.15 * ymax), 3.0);   // cap; off-scale arrowed

                // Reserved header strip + the panel grid in a TPad below, as in Step 3: the eight
                // pair-pT entries in ONE row (SetNColumns(8)) overran the canvas -- each entry's
                // text ran into the next entry's marker and the last one was clipped off the right
                // edge. Two rows of four.
                const int    kHeaderPx = 128;
                const int    canv_h    = 450 * nrow + kHeaderPx;
                const double hfrac     = (double)kHeaderPx / canv_h;
                TCanvas c(("c_s4_pteta_" + R.tag).c_str(), "", 500 * ncol, canv_h);
                auto* grid = new TPad(("grid_s4pe_" + R.tag).c_str(), "", 0., 0., 1., 1. - hfrac);
                grid->SetFillStyle(0); grid->Draw(); grid->cd(); grid->Divide(ncol, nrow);
                std::vector<TH1D*> rs_leg;   // first panel's series -> canvas-level legend
                std::vector<std::string> offscale;
                for (int iz = 1; iz <= neta; ++iz) {
                    grid->cd(iz);
                    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
                    gPad->SetTopMargin(0.10);   // label strip above the frame -- see Step 3
                    std::vector<TH1D*>& rs = rs_all[iz];
                    DrawEffFrame(0.0, R.xhi, "#DeltaR", 0.0, ymax, eps_single_text);
                    DrawUnityLine(0.0, R.xhi);
                    for (int iy = 0; iy < npt; ++iy) {
                        rs[iy]->SetMarkerStyle(ptmark[iy % ptmark.size()]);
                        rs[iy]->SetMarkerColor(ptcol[iy % ptcol.size()]);
                        rs[iy]->SetLineColor(ptcol[iy % ptcol.size()]);
                        rs[iy]->SetLineWidth(2);
                    }
                    if (iz == 1) rs_leg = rs;
                    for (int iy = 0; iy < npt; ++iy) {
                        TH1D* r = rs[iy];
                        r->Draw("E1 same");
                        for (int i = 1; i <= r->GetNbinsX(); ++i) {
                            if (r->GetBinContent(i) <= ymax) continue;
                            auto* ar = new TArrow(r->GetBinCenter(i), ymax * 0.88,
                                                  r->GetBinCenter(i), ymax * 0.985, 0.008, "|>");
                            ar->SetLineColor(ptcol[iy % ptcol.size()]);
                            ar->SetFillColor(ptcol[iy % ptcol.size()]); ar->Draw();
                            offscale.push_back(Form(
                                "#eta^{pair} #in [%.1f,%.1f), p_{T}^{pair} #in [%.1f,%.1f) GeV,"
                                " #DeltaR = %.2f: %.1f",
                                h3fn->GetZaxis()->GetBinLowEdge(iz),
                                h3fn->GetZaxis()->GetBinUpEdge(iz),
                                h3fn->GetYaxis()->GetBinLowEdge(iy + 1),
                                h3fn->GetYaxis()->GetBinUpEdge(iy + 1),
                                r->GetBinCenter(i), r->GetBinContent(i)));
                        }
                    }
                    TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.050);
                    tl.DrawLatex(0.15, 0.945, eta_label(iz).c_str());   // reserved strip, see Step 3
                }
                c.cd(0);
                auto ny = [&](double px) { return 1.0 - px / canv_h; };
                TLatex st; st.SetNDC(); st.SetTextFont(42); st.SetTextSize(22.0 / canv_h);
                st.DrawLatex(0.02, ny(26), (headline + ",  " + eps_single_text +
                             " in (p_{T}^{pair}, #eta^{pair}) cells").c_str());
                if (!rs_leg.empty()) {   // canvas-level legend, as in Step 3
                    auto* leg = new TLegend(0.02, ny(100), 0.98, ny(38));
                    leg->SetNColumns(4);
                    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(16.0 / canv_h);
                    for (int iy = 0; iy < npt; ++iy)
                        leg->AddEntry(rs_leg[iy], pt_label(iy + 1).c_str(), "lp");
                    leg->Draw();
                }
                if (!offscale.empty()) {
                    std::string note = "above the axis range: " + offscale[0];
                    for (size_t k = 1; k < offscale.size() && k < 2; ++k) note += ";  " + offscale[k];
                    if (offscale.size() > 2)
                        note += Form(";  ... (%d in total)", (int)offscale.size());
                    TLatex os_; os_.SetNDC(); os_.SetTextFont(42); os_.SetTextSize(15.0 / canv_h);
                    os_.SetTextColor(kGray + 3);
                    os_.DrawLatex(0.02, ny(118), note.c_str());
                }
                SaveCanvas(c, dir4 + "step4_eps_dr_single_" + R.tag + "_pair_eta_pt.png");
            }

            // ---- plateau + fluctuation tables (full-dR 3D, window [kPlateauLo,kPlateauHi]) ----
            using Plat = PlateauCell;   // shared with the ROOT-file writer below
            auto cell_plateau = [&](int iy, int iz) -> Plat {
                TH1D* r = cell_ratio(h3fn, h3fd, h3fA, h3fB, h3fP, h3fQ, iy, iz, neta,
                                     Form("s4plat_%s_%d_%d", sample.c_str(), iy, iz));
                auto window = [&](double lo, double hi) -> Plat {
                    double sw = 0, swv = 0; std::vector<std::pair<double,double>> vw;
                    for (int i = 1; i <= r->GetNbinsX(); ++i) {
                        const double xc = r->GetBinCenter(i);
                        if (xc < lo || xc > hi) continue;
                        const double v = r->GetBinContent(i), e = r->GetBinError(i);
                        if (e <= 0. || v == 0.) continue;
                        const double w = 1. / (e * e); sw += w; swv += w * v; vw.push_back({v, w});
                    }
                    if (sw <= 0.) return Plat{ -1, -1, -1, 0 };
                    const double mean = swv / sw, err = std::sqrt(1. / sw);
                    double swd = 0;
                    for (auto& p : vw) swd += p.second * (p.first - mean) * (p.first - mean);
                    return Plat{ mean, err, std::sqrt(swd / sw), (int)vw.size() };
                };
                Plat nom = window(kPlateauLo, kPlateauHi);
                const Plat alt = window(kPlateauSystLo, kPlateauSystHi);
                delete r;
                if (nom.nb > 0 && alt.nb > 0) nom.syst = std::fabs(nom.mean - alt.mean);
                return nom;
            };
            std::vector<std::vector<Plat>> P(npt, std::vector<Plat>(neta));
            for (int iy = 1; iy <= npt; ++iy)
                for (int iz = 1; iz <= neta; ++iz) P[iy - 1][iz - 1] = cell_plateau(iy, iz);
            auto eta_hdr = [&](int iz){ return Form("[%.1f,%.1f)",
                h3fn->GetZaxis()->GetBinLowEdge(iz), h3fn->GetZaxis()->GetBinUpEdge(iz)); };
            auto pt_hdr  = [&](int iy){ return Form("pTpair[%.1f,%.1f)",
                h3fn->GetYaxis()->GetBinLowEdge(iy), h3fn->GetYaxis()->GetBinUpEdge(iy)); };
            auto write_table_A = [&](std::ostream& os){
                os << "# Step-4 large-dR plateau eps_single (weighted mean over dR in ["
                   << kPlateauLo << "," << kPlateauHi << "]) per (pair pT, pair eta) cell.\n";
                os << "# sample=" << sample << "  WP=" << wp_text << "  " << cfg.sample_text << "\n";
                os << "# " << (is_deliverable ? "PbPb DELIVERABLE (dresses union linear terms)"
                                              : "pp VALIDATION only -- NOT applied to pp 2mu4") << "\n";
                os << "# value = plateau +- stat_error (normalized to 1; |value-1| feeds the systematic)\n";
                os << std::left << std::setw(20) << "pair-eta \\ pair-pT";
                for (int iy = 1; iy <= npt; ++iy) os << std::setw(20) << pt_hdr(iy);
                os << "\n";
                for (int iz = 1; iz <= neta; ++iz) {
                    os << std::left << std::setw(20) << eta_hdr(iz);
                    for (int iy = 1; iy <= npt; ++iy) {
                        const Plat& p = P[iy - 1][iz - 1];
                        os << std::setw(20) << (p.nb > 0 ? Form("%.4f+-%.4f", p.mean, p.err) : "--");
                    }
                    os << "\n";
                }
            };
            auto write_table_B = [&](std::ostream& os){
                os << "# Step-4 plateau FLUCTUATION per (pair pT, pair eta) cell (dR in ["
                   << kPlateauLo << "," << kPlateauHi << "]).\n";
                os << "# sample=" << sample << "  WP=" << wp_text << "\n";
                os << "# each cell = stat_err(mean) / rms_scatter / n_bins / window_syst.\n";
                os << "# window_syst = |plateau[" << kPlateauLo << "," << kPlateauHi
                   << "] - plateau[" << kPlateauSystLo << "," << kPlateauSystHi << "]|, the"
                      " plateau-window normalization systematic (-- = not evaluable).\n";
                os << std::left << std::setw(28) << "pair-eta \\ pair-pT";
                for (int iy = 1; iy <= npt; ++iy) os << std::setw(28) << pt_hdr(iy);
                os << "\n";
                for (int iz = 1; iz <= neta; ++iz) {
                    os << std::left << std::setw(28) << eta_hdr(iz);
                    for (int iy = 1; iy <= npt; ++iy) {
                        const Plat& p = P[iy - 1][iz - 1];
                        os << std::setw(28) << (p.nb > 0
                            ? Form("%.4f/%.4f/%d/%s", p.err, p.rms, p.nb,
                                   p.syst >= 0. ? Form("%.4f", p.syst) : "--")
                            : "--");
                    }
                    os << "\n";
                }
            };
            write_table_A(std::cout);
            { std::ofstream ofa(dir4 + "step4_plateau_pair_eta_pt.txt");             write_table_A(ofa); }
            { std::ofstream ofb(dir4 + "step4_plateau_fluctuation_pair_eta_pt.txt"); write_table_B(ofb); }
            std::cout << "  wrote " << dir4 << "step4_plateau_pair_eta_pt.txt + "
                      << "step4_plateau_fluctuation_pair_eta_pt.txt\n";

            // Machine-readable copy for fit_dr_corrections.cxx. UPDATE, not RECREATE: Step 3
            // created the file earlier in this same run and its maps must survive.
            WritePlateauRootFile(DrCorrPlateauFile(id, use_tight_wp), /*recreate=*/false, 4,
                                 h3fn, P, plat4.first, plat4.second, plat4_syst,
                                 sample, wp_text, "eps_single");

            // ---- ROUND 9: same-sign / opposite-sign series, exactly as Step 3 --------------
            for (const std::string& sgn : {std::string("ss"), std::string("os")}) {
                const std::string hp4 = "h_mc_single_dr_" + sgn + "_";
                TH3D* sn = GetObjOrNull<TH3D>(fmc4, hp4 + "full_vs_pt_eta_num");
                TH3D* sd = GetObjOrNull<TH3D>(fmc4, hp4 + "full_vs_pt_eta_denom");
                TH3D* sA = GetObjOrNull<TH3D>(fmc4, hp4 + "full_vs_pt_eta_errA");
                TH3D* sB = GetObjOrNull<TH3D>(fmc4, hp4 + "full_vs_pt_eta_errB");
                TH3D* sP = GetObjOrNull<TH3D>(fmc4, hp4 + "full_vs_pt_eta_covP");
                TH3D* sQ = GetObjOrNull<TH3D>(fmc4, hp4 + "full_vs_pt_eta_covQ");
                TH1D* i_n = GetObjOrNull<TH1D>(fmc4, hp4 + "full_num");
                TH1D* i_d = GetObjOrNull<TH1D>(fmc4, hp4 + "full_denom");
                TH1D* i_A = GetObjOrNull<TH1D>(fmc4, hp4 + "full_errA");
                TH1D* i_B = GetObjOrNull<TH1D>(fmc4, hp4 + "full_errB");
                TH1D* i_P = GetObjOrNull<TH1D>(fmc4, hp4 + "full_covP");
                TH1D* i_Q = GetObjOrNull<TH1D>(fmc4, hp4 + "full_covQ");
                if (!sn || !sd || !sA || !sB || !sP || !sQ ||
                    !i_n || !i_d || !i_A || !i_B || !i_P || !i_Q) {
                    std::cout << "  [Step 4] no " << hp4 << "* histograms -- sign-separated "
                              << "plateaus skipped for this sample (refill to produce them)\n";
                    continue;
                }
                auto* ri = (TH1D*)i_n->Clone(("s4incl_" + sgn).c_str());
                ri->SetDirectory(nullptr);
                ri->Divide(i_d);
                SetConditionalRatioErrors(ri, i_d, i_A, i_B, i_P, i_Q);
                const auto pi  = PlateauWeightedMean(ri, kPlateauLo, kPlateauHi);
                const auto pia = PlateauWeightedMean(ri, kPlateauSystLo, kPlateauSystHi);
                delete ri;

                std::vector<std::vector<Plat>> Ps(npt, std::vector<Plat>(neta));
                for (int iy = 1; iy <= npt; ++iy)
                    for (int iz = 1; iz <= neta; ++iz) {
                        TH1D* r = cell_ratio(sn, sd, sA, sB, sP, sQ, iy, iz, neta,
                                             Form("plat4%s_%s_%d_%d", sgn.c_str(),
                                                  sample.c_str(), iy, iz));
                        Ps[iy - 1][iz - 1] = PlateauFromRatio(r);
                        delete r;
                    }
                WritePlateauRootFile(DrCorrPlateauFile(id, use_tight_wp), /*recreate=*/false, 4,
                                     sn, Ps, pi.first, pi.second,
                                     std::fabs(pi.first - pia.first),
                                     sample, wp_text, "eps_single", sgn);
                printf("  %s (%s pairs) single-leg inclusive plateau = %.4f +- %.4f\n",
                       sample.c_str(), sgn == "ss" ? "same-sign" : "opposite-sign",
                       pi.first, pi.second);
            }
            fmc4->Close();
        }
    }

    fmc->Close(); fmc3->Close(); ffit->Close(); fdata->Close();
    std::cout << "\nplot_mc_trig_eff(" << sample << ") done.\n";
}
