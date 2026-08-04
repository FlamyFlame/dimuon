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
#include <TROOT.h>

#include <TH3D.h>

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

// =============================================================================
// CONDITIONAL (binomial-correct) ERROR ON AN INVERSE-WEIGHTED EFFICIENCY RATIO  (round 7)
//
// R(dR) = N/D with D = sum_all w  and  N = sum_fired w/eps  is an EFFICIENCY: the numerator
// is a re-weighted SUBSET of the denominator. `TH1::Divide` without option "B" propagates
//     e_R = R sqrt((e_N/N)^2 + (e_D/D)^2)
// which assumes num and den are INDEPENDENT. They are not, and the resulting bars are too
// long by ~sqrt((1+eps)/(1-eps)) -- 1.5x to 3.2x here, worst in the highest pair-pT bin where
// eps is largest. Symptom: chi2/ndf of a constant fit over the plateau ~0.25 instead of ~1.
//
// Conditioning on the MC sample (D fixed; only the Bernoulli trigger decisions fluctuate),
//     Var(N) = sum_i a_i^2 p_i (1-p_i) + 2 sum_pairs a_1 a_2 (p_12 - p_1 p_2),  a_i = w_i/eps_i
// estimated from the histograms booked in FillMCTrigEffHists.cxx as
//     Var = A - R*B  (+ covP - R^2*covQ  for Step 4, where both legs of a pair share a dR bin)
//     e_R = sqrt(Var)/D
// Unweighted limit (w=1, eps=1): A=B=N -> Var = N(1-R) -> e_R = sqrt(R(1-R)/D), as it must be.
// covP/covQ are absent for Step 3 (one entry = one pair = one Bernoulli trial) -> pass nullptr.
// =============================================================================
void SetConditionalRatioErrors(TH1D* r, const TH1D* den, const TH1D* A, const TH1D* B,
                               const TH1D* covP = nullptr, const TH1D* covQ = nullptr)
{
    for (int i = 1; i <= r->GetNbinsX(); ++i) {
        const double D = den->GetBinContent(i);
        if (D <= 0.) { r->SetBinError(i, 0.); continue; }
        const double R  = r->GetBinContent(i);
        const double a  = A->GetBinContent(i), b = B->GetBinContent(i);
        const double cp = covP ? covP->GetBinContent(i) : 0.;
        const double cq = covQ ? covQ->GetBinContent(i) : 0.;
        double var = a - R * b + cp - R * R * cq;

        // BOUNDARY CASE. When every effective entry in the bin fired (p_i -> 1) the conditional
        // binomial variance genuinely vanishes -- the k=n binomial artefact -- and `var` comes
        // out as 0, slightly negative, or a catastrophic cancellation of terms many orders of
        // magnitude larger (seen in the near-empty high-pair-pT cells of the 10k-event overlay
        // TEST sample: A=1.1e-06 vs var=1e-22). A ~0 error is NOT safe: the plateau weighted
        // mean weights by 1/e^2, so such a bin would either be dropped (e=0) or completely
        // dominate the mean (e=1e-7). Detect it by comparing var with the SCALE of the terms
        // that built it, and fall back to the "1/n rule" (the 68% bound for k=n) with the
        // effective denominator count n_eff = (D/e_D)^2, e_D = sqrt(sum w^2) from Sumw2.
        const double scale = a + R * b + std::fabs(cp) + R * R * cq;
        if (var > 1e-6 * scale) { r->SetBinError(i, std::sqrt(var) / D); continue; }

        const double eD = den->GetBinError(i);
        const double neff = (eD > 0.) ? (D / eD) * (D / eD) : 1.0;
        r->SetBinError(i, (neff > 0.) ? R / neff : 0.);
    }
}

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

    SampleCfg c;
    if (sample == "pp") {
        c.mc_dir      = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/";
        c.mc_label    = "pp24";
        c.data_file   = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
                        "histograms_real_pairs_pp_2024_single_mu4_fine_q_eta_bin"
                        + data_wp + ".root";
        c.ctr         = "";
        c.out_base    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                        "pp_trigger_efficiency/mc_based/";
        c.sample_text = "Pythia8 pp24 fullsim";
        c.data_text   = "pp24 data";
        c.eps_dr_text = "#varepsilon_{#DeltaR}^{2mu4}";
        c.step2_coarse = false;  // pp has ~4x the pair statistics: the fine axes are readable
    } else if (sample == "pp_full") {
        // pp24 FULL sample. Physics identical to "pp" -- SAME pp24 data reference, SAME
        // 2mu4 product weighting, SAME plot LOCATION (this is the canonical pp trig-eff
        // deliverable, which the full sample now supersedes -- back up the TEST-sample plots
        // first, done by the pipeline / the run wrapper). ONLY the MC inputs differ:
        // mc_dir = full-sample dir, mc_label = "pp24_full" (reads the _full intermediate hists,
        // so the hists/fits never clobber the TEST ones).
        c.mc_dir      = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/";
        c.mc_label    = "pp24_full";
        c.data_file   = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
                        "histograms_real_pairs_pp_2024_single_mu4_fine_q_eta_bin"
                        + data_wp + ".root";
        c.ctr         = "";
        c.out_base    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                        "pp_trigger_efficiency/mc_based/";
        c.sample_text = "Pythia8 pp24 fullsim (FULL sample)";
        c.data_text   = "pp24 data";
        c.eps_dr_text = "#varepsilon_{#DeltaR}^{2mu4}";
        c.step2_coarse = false;  // the full sample has far MORE pair statistics than the test
    } else if (sample == "overlay") {
        c.mc_dir      = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/";
        c.mc_label    = "hijing_overlay_pbpb23";
        c.data_file   = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/"
                        "histograms_real_pairs_pbpb_2023_single_mu4_fine_q_eta_bin"
                        + data_wp + ".root";
        c.ctr         = "_ctr0_5";   // D2: overlay compares ONLY to PbPb23 data 0-5%
        c.out_base    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                        "pbpb_trigger_efficiency/mc_based/";
        c.sample_text = "HIJING overlay Pb+Pb23 cond., 0-5%";
        c.data_text   = "Pb+Pb23 data 0-5%";
        c.eps_dr_text = "#varepsilon_{#DeltaR}^{cross}";
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
        c.mc_dir      = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_no_overlay_test_sample/";
        c.mc_label    = "r17663_no_overlay";
        c.data_file   = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
                        "histograms_real_pairs_pp_2024_single_mu4_fine_q_eta_bin"
                        + data_wp + ".root";
        c.ctr         = "";          // no centrality: there is no overlaid event
        c.out_base    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                        "r17663_no_overlay_trigger_efficiency/mc_based/";
        c.sample_text = "Pythia8 pp, r17663 (HI cond., no overlay)";
        c.data_text   = "pp24 data";
        c.eps_dr_text = "#varepsilon_{#DeltaR}^{2mu4}";
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
        c.cmp_text     = "MC Pb+Pb23 cond., WITH HIJING overlay (r17618)";
        // Black series on the q.eta panels = pp24-CONDITIONS fullsim MC (FULL sample), NOT
        // pp data: the r17663 study is an MC-vs-MC comparison of reco-tag CONFIGURATIONS, so
        // all three curves are MC (red = Pb+Pb23 cond. no overlay; black = pp24 cond.;
        // blue = Pb+Pb23 cond. with HIJING overlay). Full sample chosen for statistics.
        // FitMCSinglesEffcy writes the unqualified basename, distinguished here only by the
        // full_sample directory (Remaining Work 3c).
        c.qeta_black_mc_file = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/"
                               "single_mu_effcy_pT_fit_mc"
                             + std::string(use_tight_wp ? "" : "_medium_wp") + ".root";
        c.qeta_black_mc_text = "MC pp24 cond. (fullsim, full sample)";
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

// fine q.eta bins (CommonEffcyConfig.h q_eta_proj_ranges_fine_excl_gap). The forward
// (-2.4,-2.0) bin is split into (-2.4,-2.2)+(-2.2,-2.0) (round-5 change #2) -> 11 bins.
const std::vector<std::string> kQEtaSuffix = {
    "minus2_40_TO_minus2_20", "minus2_20_TO_minus2_00",
    "minus2_00_TO_minus1_60", "minus1_60_TO_minus1_30",
    "minus0_90_TO_minus0_50", "minus0_50_TO_minus0_10", "0_10_TO_0_50",
    "0_50_TO_1_00", "1_30_TO_1_60", "1_60_TO_2_00", "2_00_TO_2_20"};
const std::vector<std::pair<double,double>> kQEtaRange = {
    {-2.4,-2.2},{-2.2,-2.0},{-2.0,-1.6},{-1.6,-1.3},{-0.9,-0.5},{-0.5,-0.1},
    { 0.1, 0.5},{ 0.5, 1.0},{ 1.3, 1.6},{ 1.6, 2.0},{ 2.0, 2.2}};

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

// Step-3 plateau window: well-separated muons must decorrelate, so eps_dR is flat here
// (§3.3 diagnostic 1). [1,4] rather than [1,3] -- the extra separation is still plateau,
// and the overlay needs every pair it can get.
const double kPlateauLo = 1.0;
const double kPlateauHi = 4.0;

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

} // namespace

// ================================================================= main

void plot_mc_trig_eff(const std::string& sample = "pp", bool use_tight_wp = true)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gErrorIgnoreLevel = kWarning;

    const SampleCfg cfg = MakeCfg(sample, use_tight_wp);

    // WP config (registry: Analysis/docs/muon_wp_registry.md): Tight nominal unsuffixed;
    // Medium inputs carry _medium_wp. BOTH the MC inputs AND the data tag-and-probe file are
    // WP-keyed (see MakeCfg) -- the data Medium reference is produced by running the data RDF
    // with isTight=false. Never compare across working points.
    const std::string wp_suf  = use_tight_wp ? "" : "_medium_wp";
    const std::string wp_text = use_tight_wp ? "Tight muons" : "Medium muons";
    const std::string headline = cfg.sample_text + ", " + wp_text;

    TFile* fmc   = OpenFile(cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf + ".root");
    TFile* fmc3  = OpenFile(cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf + "_step3.root");
    TFile* ffit  = OpenFile(cfg.mc_dir + "single_mu_effcy_pT_fit_mc" + wp_suf + ".root");
    TFile* fdata = OpenFile(cfg.data_file);

    // Medium-WP plots go into a medium/ SUBDIRECTORY of each step dir (same filenames as
    // Tight); Tight (nominal) stays in the step dir root (user request, 2026-07-16).
    const std::string wp_dir = use_tight_wp ? "" : "medium/";
    const std::string dir1 = cfg.out_base + "step1_singles_data_mc/"  + wp_dir;
    // Step-2 is produced 3x (round-5 #3): the L1, HLT|L1 and full-chain efficiencies each go
    // into their OWN subdirectory of step2_dr_binned_singles/ (built inside the stage loop).
    const std::string dir2_base = cfg.out_base + "step2_dr_binned_singles/";
    const std::string dir3 = cfg.out_base + "step3_dr_correction/"     + wp_dir;
    for (const auto& d : {dir1, dir3}) gSystem->mkdir(d.c_str(), kTRUE);

    const std::string mc_leg   = "MC direct P(mu4 | reco #mu)";
    const std::string data_leg = cfg.data_text + " T&P P(2mu4 | mu4 tag, #DeltaR>0.8)";

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
                           qeta_all_mc ? "r17663 / pp24" : "MC / data", 0.5, 2.6);
            auto* grat = DivideGraphClean(gmc, gda);
            StyleGraph(grat, kMCColor, 21, 0.7);
            grat->Draw("PZ same");
            MarkOffScale(grat, 0.5, 2.6, kMCColor);
            c.cd(static_cast<int>(iq) + 1);
        }
        // legend / label pad
        c.cd(12);  // legend/label pad: 11 q.eta bins now fill pads 1-11 (round-5 #2)
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
            leg->SetTextSize(0.037);
            leg->SetMargin(0.12);
            leg->AddEntry(gm, "MC Pb+Pb23 cond., no HIJING overlay (+fit)", "lp");
            leg->AddEntry(gd, "MC pp24 conditions (fullsim, full sample)", "lp");
            if (fcmp) {
                auto* gc = new TGraphAsymmErrors(); StyleGraph(gc, kCmpColor, 22);
                leg->AddEntry(gc, "MC Pb+Pb23 cond., with HIJING overlay", "lp");
            }
            leg->Draw();
            TLatex note;
            note.SetNDC();
            note.SetTextSize(0.036);
            note.SetTextFont(42);
            note.DrawLatex(0.02, 0.42, "All three: MC direct P(mu4 | reco #mu), no tag-and-probe.");
            note.DrawLatex(0.02, 0.34, "red = r17663 (no overlay), blue = r17618 (with overlay):");
            note.DrawLatex(0.02, 0.28, "SAME reco tag family, differ ONLY by the HIJING overlay.");
            note.DrawLatex(0.02, 0.18, "Ratio pad: r17663 / pp24-cond. MC (both #mu^{#pm} summed");
            note.DrawLatex(0.02, 0.12, "into the panel's q#upoint#eta bin).");
        } else {
            auto* leg = new TLegend(0.02, 0.45, 0.98, 0.80);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.05);
            leg->AddEntry(gm, (mc_leg + " (+ fit)").c_str(), "lp");
            leg->AddEntry(gd, (cfg.data_text + " tag&probe").c_str(), "lp");
            if (fcmp) {
                auto* gc = new TGraphAsymmErrors(); StyleGraph(gc, kCmpColor, 22);
                leg->AddEntry(gc, (cfg.cmp_text + ", no fit").c_str(), "lp");
            }
            leg->Draw();
            TLatex note;
            note.SetNDC();
            note.SetTextSize(0.045);
            note.SetTextFont(42);
            note.DrawLatex(0.02, 0.30, "Data: T&P P(2mu4 | mu4 tag, #DeltaR>0.8 pairs);");
            note.DrawLatex(0.02, 0.22, "MC: direct conditional, no T&P");
        }
        SaveCanvas(c, dir1 + "step1_eff_pt_in_q_eta_bins_" + kCharges[ic] + ".png");
    }

    // ================================================================
    // Step 2 (§3.2): DeltaR-binned MC singles efficiency
    // ================================================================
    std::cout << "\n===== Step 2 (" << sample << ", " << wp_text << ") =====\n";

    const std::string leg_incl = "inclusive singles (Step 1)";

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
    struct Step2Stage { std::string sub, numk, denk, eff; };
    const std::vector<Step2Stage> step2_stages = {
        {"full_chain", "num",    "denom", "P(full mu4 chain | reco #mu)"},
        {"L1",         "numl1",  "denom", "P(L1 MU3V | reco #mu)"},
        {"HLT",        "numhlt", "numl1", "P(HLT | fires L1, reco #mu)"}};
    for (const auto& st : step2_stages) {
      const std::string dir2 = dir2_base + st.sub + "/" + wp_dir;
      gSystem->mkdir(dir2.c_str(), kTRUE);
      const std::string NUMK = st.numk, DENK = st.denk;
      const std::string stage_tag = " [" + st.sub + "]";
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
            leg->SetFillColorAlpha(kWhite, 0.75);
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
            DrawHeadline(headline + ", " + kChargeTex[ic] + stage_tag, 0.14, 0.955, 0.05);

            // R3 ratio pad: §3.2 asks whether the three DeltaR series AGREE at fixed
            // kinematics. Ratio to the isolated (DeltaR >= 1.0) series = that test.
            pads.second->cd();
            DrawRatioFrame(v.xlo, v.xhi, v.xtitle, "/ #DeltaR #geq 1", 0.4, 1.9);
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
            DrawRatioFrame(4.0, 60.0, "p_{T} [GeV]", "/ #DeltaR #geq 1", 0.4, 1.9);
            for (size_t id = 0; id + 1 < gdr.size(); ++id) {
                auto* g = DivideGraphClean(gdr[id], gdr.back());
                StyleGraph(g, kDrColor[id], kDrMarker[id], 0.7);
                g->Draw("PZ same");
                MarkOffScale(g, 0.4, 1.9, kDrColor[id]);
            }
            c.cd(static_cast<int>(iq) + 1);
        }
        c.cd(12);  // legend/label pad: 11 q.eta bins now fill pads 1-11 (round-5 #2)
        DrawHeadline(headline + ", " + kChargeTex[ic] + stage_tag, 0.02, 0.88, 0.048); // 0.048: long overlay headline + charge must fit (review iter 1)
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

            struct Var { std::string key, tex; Color_t col; Style_t mk; };
            std::vector<Var> vars = {
                {"orig", "original MC (round-7 selection)",       kBlack,     20},
                {"vtx",  "+ 1 reconstructed vertex",              kBlue + 1,  21},
                {"ptm",  "+ |#Deltap_{T}|/p_{T}^{truth} < thr",   kGreen + 2, 22},
                {"both", "+ both",                                kRed + 1,   23}};

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
                    else { dropped.push_back(v.key);
                           std::cout << "  [sanity] variant '" << v.key
                                     << "' has an EMPTY denominator -> requirement inapplicable "
                                        "to this sample; series dropped\n"; }
                }
                vars = keep;
            }
            const std::string drop_note = dropped.empty() ? std::string()
                : ("inapplicable to this sample (empty): " +
                   [&]{ std::string s; for (size_t i = 0; i < dropped.size(); ++i)
                        s += (i ? ", " : "") + dropped[i]; return s; }());

            auto eff_of = [&](const std::string& base, const std::string& chg,
                              const std::string& v) -> TGraphAsymmErrors* {
                TH1D* n = GetObj<TH1D>(fsan, "h_sanity_" + base + "_num_"   + chg + "_" + v);
                TH1D* d = GetObj<TH1D>(fsan, "h_sanity_" + base + "_denom_" + chg + "_" + v);
                return BayesEff(n, d);
            };

            // ---- (a) eff vs pT and vs q.eta, one pad per charge, ratio pad vs `orig` ----
            struct Obs { std::string base, xt; double xlo, xhi; bool logx; };
            const std::vector<Obs> obs = {{"pt", "p_{T} [GeV]", 4.0, 60.0, true},
                                          {"q_eta", "q#eta", -2.4, 2.4, false}};
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
                    DrawEffFrame(O.xlo, O.xhi, O.xt, 0.90, 1.10, "variant / original");
                    DrawUnityLine(O.xlo, O.xhi);
                    for (size_t i = 1; i < gs.size(); ++i) {   // skip `orig` (ratio 1 by construction)
                        auto* gr = DivideGraphClean(gs[i], gs[0]);
                        StyleGraph(gr, vars[i].col, vars[i].mk, 0.8);
                        gr->Draw("PZ same");
                    }
                }
                DrawHeadline(headline + "  --  Step-1 sanity check");
                if (!drop_note.empty()) {
                    c.cd(0);
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
                    tl.DrawLatex(0.20, 0.90, Form("%.2f < q#eta < %.2f",
                                                  kQEtaRange[iq].first, kQEtaRange[iq].second));
                    if (iq == 0) {
                        auto* leg = new TLegend(0.30, 0.14, 0.95, 0.40);
                        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
                        for (size_t i = 0; i < vars.size(); ++i)
                            leg->AddEntry(gs[i], vars[i].tex.c_str(), "lp");
                        leg->Draw();
                    }
                }
                DrawHeadline(headline + "  --  Step-1 sanity, " + kChargeTex[ic]);
                if (!drop_note.empty()) {
                    c.cd(0);
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
    printf("  %s plateau (weighted mean, dR in [%.0f,%.0f]): %.4f +- %.4f\n",
           sample.c_str(), kPlateauLo, kPlateauHi, plateau.first, plateau.second);

    auto DrawStep3 = [&](TH1D* r, double xlo, double xhi, const std::string& png,
                         bool plateau_in_range) {
        TCanvas c(("c_" + png).c_str(), "", 900, 700);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        double ymax = 0.;
        for (int i = 1; i <= r->GetNbinsX(); ++i)
            ymax = std::max(ymax, r->GetBinContent(i) + r->GetBinError(i));
        ymax = std::max(1.15, 1.15 * ymax);   // y linear, auto but include 1
        DrawEffFrame(xlo, xhi, "#DeltaR", 0.0, ymax, cfg.eps_dr_text);
        DrawUnityLine(xlo, xhi);

        // fitted plateau line: solid over the plateau window [kPlateauLo,kPlateauHi] where in range,
        // dotted across the pad otherwise (zoom canvas)
        if (plateau_in_range) {
            auto* lp = new TLine(kPlateauLo, plateau.first, kPlateauHi, plateau.first);
            lp->SetLineColor(kBlue + 1);
            lp->SetLineWidth(3);
            lp->Draw("same");
        } else {
            auto* lp = new TLine(xlo, plateau.first, xhi, plateau.first);
            lp->SetLineColor(kBlue + 1);
            lp->SetLineWidth(2);
            lp->SetLineStyle(3);
            lp->Draw("same");
        }

        r->SetMarkerStyle(20);
        r->SetMarkerColor(kMCColor);
        r->SetLineColor(kMCColor);
        r->SetLineWidth(2);
        r->Draw("E1 same");

        DrawHeadline(headline);
        TLatex tl;
        tl.SetNDC();
        tl.SetTextFont(42);
        tl.SetTextSize(0.035);
        tl.DrawLatex(0.40, 0.86, Form("plateau #LT#DeltaR#in[%.0f,%.0f]#GT = %.3f #pm %.3f",
                                      kPlateauLo, kPlateauHi, plateau.first, plateau.second));
        tl.DrawLatex(0.40, 0.80, (cfg.eps_dr_text +
            " = P(trig | #DeltaR) / (#varepsilon_{1}#varepsilon_{2}), MC #varepsilon in weights").c_str());
        SaveCanvas(c, dir3 + png + ".png");
    };
    DrawStep3(r_zoom, 0.0, 1.0,  "step3_eps_dr_zoom", false);
    DrawStep3(r_full, 0.0, 5.75, "step3_eps_dr_full", true);

    // --- pair-pT-binned zoom ratio (4 slices of the pT_bins_120 axis) ---
    {
        TH2D* h2n = GetObj<TH2D>(fmc3, "h_mc_dr_zoom_vs_pair_pt_num");
        TH2D* h2d = GetObj<TH2D>(fmc3, "h_mc_dr_zoom_vs_pair_pt_denom");
        TH2D* h2A = GetObj<TH2D>(fmc3, "h_mc_dr_zoom_vs_pair_pt_errA");
        TH2D* h2B = GetObj<TH2D>(fmc3, "h_mc_dr_zoom_vs_pair_pt_errB");
        // slice edges aligned to the pair-pT axis bin edges (F2)
        const std::vector<std::pair<int,int>> ybins = {{1,3},{4,6},{7,9},{10,15}};
        // last slice: bright kMagenta, NOT kMagenta+2 -- the darkened shade reads as another
        // dark red/blue against kRed+1 / kBlue+1 and the series cannot be told apart
        const std::vector<Color_t> scol   = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta};
        const std::vector<Style_t> smark  = {20, 21, 22, 23};

        TCanvas c("c_step3_ptslices", "", 900, 700);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);

        std::vector<TH1D*> ratios;
        std::vector<std::string> slabels;
        double ymax = 0.;
        for (size_t is = 0; is < ybins.size(); ++is) {
            TH1D* n = h2n->ProjectionX(Form("s3_n_%zu", is), ybins[is].first, ybins[is].second);
            TH1D* d = h2d->ProjectionX(Form("s3_d_%zu", is), ybins[is].first, ybins[is].second);
            TH1D* a = h2A->ProjectionX(Form("s3_a_%zu", is), ybins[is].first, ybins[is].second);
            TH1D* b = h2B->ProjectionX(Form("s3_b_%zu", is), ybins[is].first, ybins[is].second);
            auto* r = (TH1D*)n->Clone(Form("s3_r_%zu", is));
            r->SetDirectory(nullptr);
            r->Divide(d);
            SetConditionalRatioErrors(r, d, a, b);
            delete a; delete b;
            ratios.push_back(r);
            const double plo = h2n->GetYaxis()->GetBinLowEdge(ybins[is].first);
            const double phi = h2n->GetYaxis()->GetBinUpEdge(ybins[is].second);
            slabels.push_back(Form("%.1f < p_{T}^{pair} < %.1f GeV", plo, phi));
            for (int i = 1; i <= r->GetNbinsX(); ++i)
                ymax = std::max(ymax, r->GetBinContent(i) + r->GetBinError(i));
        }
        // CAP the auto-range: max+error is set by the noisiest high-pair-pT bin (the overlay
        // has one at 6.4 +- 5.5), which would squeeze all four series -- and the small-DeltaR
        // structure this panel exists to show -- into the bottom sliver of the pad.
        // A cap hides points, so every point pushed off-scale is MARKED with an up-arrow and
        // listed on the canvas: silently dropping data from a physics figure is not allowed.
        ymax = std::min(std::max(1.15, 1.15 * ymax), 3.0);
        DrawEffFrame(0.0, 1.0, "#DeltaR", 0.0, ymax, cfg.eps_dr_text);
        DrawUnityLine(0.0, 1.0);

        auto* leg = new TLegend(0.45, 0.68, 0.92, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        std::vector<std::string> offscale;
        for (size_t is = 0; is < ratios.size(); ++is) {
            ratios[is]->SetMarkerStyle(smark[is]);
            ratios[is]->SetMarkerColor(scol[is]);
            ratios[is]->SetLineColor(scol[is]);
            ratios[is]->SetLineWidth(2);
            ratios[is]->Draw("E1 same");
            leg->AddEntry(ratios[is], slabels[is].c_str(), "lp");

            // mark every point whose CENTRAL VALUE is above the capped frame
            for (int i = 1; i <= ratios[is]->GetNbinsX(); ++i) {
                const double v = ratios[is]->GetBinContent(i);
                if (v <= ymax) continue;
                const double x = ratios[is]->GetBinCenter(i);
                auto* ar = new TArrow(x, ymax * 0.88, x, ymax * 0.985, 0.012, "|>");
                ar->SetLineColor(scol[is]);
                ar->SetFillColor(scol[is]);
                ar->SetLineWidth(2);
                ar->Draw();
                offscale.push_back(Form("#DeltaR=%.2f: %.1f #pm %.1f", x, v,
                                        ratios[is]->GetBinError(i)));
            }
        }
        leg->Draw();
        DrawHeadline(headline);
        if (!offscale.empty()) {
            // BELOW the legend (which occupies y 0.68-0.88): the note used to be drawn at
            // y=0.86, straight through the legend box, making both unreadable -- and this
            // note is precisely what keeps the y-cap honest. Wrapped 2 entries per line.
            TLatex note;
            note.SetNDC();
            note.SetTextFont(42);
            note.SetTextSize(0.026);
            note.SetTextColor(kGray + 3);
            double y = 0.63;
            for (size_t i = 0; i < offscale.size(); i += 2) {
                std::string txt = (i == 0) ? "above scale (arrows): " : "  ";
                txt += offscale[i];
                if (i + 1 < offscale.size()) txt += ", " + offscale[i + 1];
                note.DrawLatex(0.45, y, txt.c_str());
                y -= 0.035;
            }
        }
        SaveCanvas(c, dir3 + "step3_eps_dr_zoom_pair_pt_slices.png");
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

        auto pt_label  = [&](int iy){ return std::string(Form("%.0f < p_{T}^{pair} < %.0f GeV",
            h3fn->GetYaxis()->GetBinLowEdge(iy), h3fn->GetYaxis()->GetBinUpEdge(iy))); };
        auto eta_label = [&](int iz){ return std::string(Form("%.1f < #eta^{pair} < %.1f",
            h3fn->GetZaxis()->GetBinLowEdge(iz), h3fn->GetZaxis()->GetBinUpEdge(iz))); };

        // eps_dR(dR) for one (pt bin iy, eta bin iz) cell. Central value = num/denom;
        // ERROR = conditional/binomial-correct form (SetConditionalRatioErrors).
        auto cell_ratio = [](TH3D* hn, TH3D* hd, TH3D* ha, TH3D* hb,
                             int iy, int iz, const char* nm) -> TH1D* {
            TH1D* n = hn->ProjectionX(Form("%s_n", nm), iy, iy, iz, iz, "e");
            TH1D* d = hd->ProjectionX(Form("%s_d", nm), iy, iy, iz, iz, "e");
            TH1D* a = ha->ProjectionX(Form("%s_a", nm), iy, iy, iz, iz, "e");
            TH1D* b = hb->ProjectionX(Form("%s_b", nm), iy, iy, iz, iz, "e");
            auto* r = (TH1D*)n->Clone(nm);
            r->SetDirectory(nullptr);
            r->Divide(d);
            SetConditionalRatioErrors(r, d, a, b);
            delete n; delete d; delete a; delete b;
            return r;
        };

        const std::vector<Color_t> ptcol  = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta};
        const std::vector<Style_t> ptmark = {20, 21, 22, 23};

        // ---- two panel plots (zoom + full dR) ----
        struct Rng { TH3D* n; TH3D* d; TH3D* a; TH3D* b; std::string tag; double xhi; };
        const std::vector<Rng> rngs = {{h3zn, h3zd, h3zA, h3zB, "zoom", 1.0},
                                       {h3fn, h3fd, h3fA, h3fB, "full", 5.75}};
        // subplot grid: nrows >= ncols, nrows ~ sqrt(neta) (feedback_subplot_layout)
        const int ncol = (int)std::ceil(std::sqrt((double)neta));
        const int nrow = (int)std::ceil((double)neta / ncol);
        for (const auto& R : rngs) {
            TCanvas c(("c_s3_pteta_" + R.tag).c_str(), "", 500 * ncol, 450 * nrow);
            c.Divide(ncol, nrow);
            for (int iz = 1; iz <= neta; ++iz) {
                c.cd(iz);
                gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
                std::vector<TH1D*> rs;
                double ymax = 0.;
                for (int iy = 1; iy <= npt; ++iy) {
                    TH1D* r = cell_ratio(R.n, R.d, R.a, R.b, iy, iz,
                                         Form("s3pe_%s_%s_%d_%d", sample.c_str(), R.tag.c_str(), iy, iz));
                    rs.push_back(r);
                    for (int i = 1; i <= r->GetNbinsX(); ++i)
                        ymax = std::max(ymax, r->GetBinContent(i) + r->GetBinError(i));
                }
                ymax = std::min(std::max(1.15, 1.15 * ymax), 3.0);   // cap; off-scale points arrowed
                DrawEffFrame(0.0, R.xhi, "#DeltaR", 0.0, ymax, cfg.eps_dr_text);
                DrawUnityLine(0.0, R.xhi);
                for (int iy = 0; iy < npt; ++iy) {
                    TH1D* r = rs[iy];
                    r->SetMarkerStyle(ptmark[iy % ptmark.size()]);
                    r->SetMarkerColor(ptcol[iy % ptcol.size()]);
                    r->SetLineColor(ptcol[iy % ptcol.size()]);
                    r->SetLineWidth(2);
                    r->Draw("E1 same");
                    for (int i = 1; i <= r->GetNbinsX(); ++i) {
                        if (r->GetBinContent(i) <= ymax) continue;
                        auto* ar = new TArrow(r->GetBinCenter(i), ymax * 0.88,
                                              r->GetBinCenter(i), ymax * 0.985, 0.008, "|>");
                        ar->SetLineColor(ptcol[iy % ptcol.size()]);
                        ar->SetFillColor(ptcol[iy % ptcol.size()]);
                        ar->Draw();
                    }
                }
                // eta-bin label INSIDE the frame top-left (was in the top margin at y=0.94,
                // where it collided with the canvas super-title on the top pad row -- review WARNING).
                TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.050);
                tl.DrawLatex(0.17, 0.86, eta_label(iz).c_str());
                if (iz == 1) {   // pair-pT-bin legend, once
                    auto* leg = new TLegend(0.42, 0.66, 0.92, 0.90);
                    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
                    for (int iy = 0; iy < npt; ++iy) leg->AddEntry(rs[iy], pt_label(iy + 1).c_str(), "lp");
                    leg->Draw();
                }
            }
            c.cd(0);
            // Super-title in the empty top strip. Use a RAW TLatex (NOT DrawHeadline, which floors
            // the size at 0.030 -> the long title then overflows the right edge; review WARNING).
            // At 0.016 the longest title (r17663 medium/zoom) fits the 1500px width, and baseline
            // 0.978 + ~0.016 height < 1.0 clears the top edge; it sits above the in-frame eta labels.
            TLatex st; st.SetNDC(); st.SetTextFont(42); st.SetTextSize(0.016);
            st.DrawLatex(0.03, 0.978, (headline + "  --  " + R.tag + " #DeltaR, " + cfg.eps_dr_text +
                         " in (p_{T}^{pair}, #eta^{pair}) cells").c_str());
            SaveCanvas(c, dir3 + "step3_eps_dr_" + R.tag + "_pair_eta_pt.png");
        }

        // ---- tables (full-dR 3D, plateau window [kPlateauLo, kPlateauHi]) ----
        struct Plat { double mean, err, rms; int nb; };
        auto cell_plateau = [&](int iy, int iz) -> Plat {
            TH1D* r = cell_ratio(h3fn, h3fd, h3fA, h3fB, iy, iz,
                                 Form("plat_%s_%d_%d", sample.c_str(), iy, iz));
            double sw = 0, swv = 0;
            std::vector<std::pair<double,double>> vw;  // (value, weight)
            for (int i = 1; i <= r->GetNbinsX(); ++i) {
                const double xc = r->GetBinCenter(i);
                if (xc < kPlateauLo || xc > kPlateauHi) continue;
                const double v = r->GetBinContent(i), e = r->GetBinError(i);
                if (e <= 0. || v == 0.) continue;
                const double w = 1. / (e * e); sw += w; swv += w * v; vw.push_back({v, w});
            }
            delete r;
            if (sw <= 0.) return Plat{ -1, -1, -1, 0 };
            const double mean = swv / sw, err = std::sqrt(1. / sw);
            double swd = 0; for (auto& p : vw) swd += p.second * (p.first - mean) * (p.first - mean);
            const double rms = std::sqrt(swd / sw);      // weighted RMS scatter about the mean
            return Plat{ mean, err, rms, (int)vw.size() };
        };

        std::vector<std::vector<Plat>> P(npt, std::vector<Plat>(neta));
        for (int iy = 1; iy <= npt; ++iy)
            for (int iz = 1; iz <= neta; ++iz) P[iy - 1][iz - 1] = cell_plateau(iy, iz);

        // Table A: plateau value +- stat error ; Table B: fluctuation (stat err + RMS scatter)
        auto eta_hdr = [&](int iz){ return Form("[%.1f,%.1f)",
            h3fn->GetZaxis()->GetBinLowEdge(iz), h3fn->GetZaxis()->GetBinUpEdge(iz)); };
        auto pt_hdr  = [&](int iy){ return Form("pTpair[%.0f,%.0f)",
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
                  " bins about their mean) / n_bins.\n";
            os << "# rms_scatter measures how flat/noisy the plateau is; compare |mean-1| to it +"
                  " stat_err to judge whether the deviation is significant (systematic size).\n";
            os << std::left << std::setw(24) << "pair-eta \\ pair-pT";
            for (int iy = 1; iy <= npt; ++iy) os << std::setw(24) << pt_hdr(iy);
            os << "\n";
            for (int iz = 1; iz <= neta; ++iz) {
                os << std::left << std::setw(24) << eta_hdr(iz);
                for (int iy = 1; iy <= npt; ++iy) {
                    const Plat& p = P[iy - 1][iz - 1];
                    os << std::setw(24) << (p.nb > 0 ? Form("%.4f/%.4f/%d", p.err, p.rms, p.nb) : "--");
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
    }

    // ================================================================
    // Step 4 (§3.4): single-leg ΔR correction ε_single(ΔR) = inverse-weighted num / denom.
    //   Leg-level analog of Step 3. Deliverable for the PbPb UNION linear terms (overlay);
    //   for pp it is a VALIDATION only (2mu4 product absorbs the single-leg ΔR into ε_ΔR^2mu4)
    //   -- annotated on the pp panels. Graceful skip if _step4.root is absent.
    // ================================================================
    {
        const std::string step4_path =
            cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf + "_step4.root";
        TFile* fmc4 = TFile::Open(step4_path.c_str(), "READ");
        if (!fmc4 || fmc4->IsZombie()) {
            std::cout << "\n[Step 4] no " << step4_path << " -- skipping single-leg ΔR plots "
                      << "(run FillMCTrigEffHists do_step4=true first).\n";
        } else {
            std::cout << "\n===== Step 4 (" << sample << ", " << wp_text << ") =====\n";
            const std::string dir4 = cfg.out_base + "step4_dr_correction_singles/" + wp_dir;
            gSystem->mkdir(dir4.c_str(), kTRUE);
            const std::string eps_single_text = "#varepsilon_{#DeltaR}^{single}";
            const bool is_deliverable = (sample == "overlay");   // PbPb union; else validation
            const std::string val_note = is_deliverable
                ? std::string("dresses the union linear terms (#varepsilon_{1}+#varepsilon_{2})")
                : std::string("VALIDATION only -- NOT applied to pp 2mu4");

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
            printf("  %s single-leg plateau (weighted mean, dR in [%.0f,%.0f]): %.4f +- %.4f\n",
                   sample.c_str(), kPlateauLo, kPlateauHi, plat4.first, plat4.second);

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
                if (plateau_in_range) {
                    auto* lp = new TLine(kPlateauLo, plat4.first, kPlateauHi, plat4.first);
                    lp->SetLineColor(kBlue + 1); lp->SetLineWidth(3); lp->Draw("same");
                } else {
                    auto* lp = new TLine(xlo, plat4.first, xhi, plat4.first);
                    lp->SetLineColor(kBlue + 1); lp->SetLineWidth(2); lp->SetLineStyle(3);
                    lp->Draw("same");
                }
                r->SetMarkerStyle(20);
                r->SetMarkerColor(kMCColor);
                r->SetLineColor(kMCColor);
                r->SetLineWidth(2);
                r->Draw("E1 same");
                DrawHeadline(headline);
                TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.035);
                tl.DrawLatex(0.40, 0.86, Form("plateau #LT#DeltaR#in[%.0f,%.0f]#GT = %.3f #pm %.3f",
                                              kPlateauLo, kPlateauHi, plat4.first, plat4.second));
                tl.DrawLatex(0.40, 0.80, (eps_single_text +
                    " = P(#mu fires | #DeltaR) / #varepsilon(p_{T},q#eta), MC #varepsilon in weights").c_str());
                tl.SetTextSize(0.030);
                tl.SetTextColor(is_deliverable ? (kGray + 3) : (kRed + 2));
                tl.DrawLatex(0.40, 0.74, val_note.c_str());
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
            auto pt_label  = [&](int iy){ return std::string(Form("%.0f < p_{T}^{pair} < %.0f GeV",
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
            const std::vector<Color_t> ptcol  = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta};
            const std::vector<Style_t> ptmark = {20, 21, 22, 23};

            // pair-pT slices (integrate eta): zoom dR
            {
                TCanvas c("c_step4_ptslices", "", 900, 700);
                gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
                std::vector<TH1D*> rs; double ymax = 0.;
                for (int iy = 1; iy <= npt; ++iy) {
                    TH1D* r = cell_ratio(h3zn, h3zd, h3zA, h3zB, h3zP, h3zQ, iy, 0, neta,
                                         Form("s4pt_%s_%d", sample.c_str(), iy));
                    rs.push_back(r);
                    for (int i = 1; i <= r->GetNbinsX(); ++i)
                        ymax = std::max(ymax, r->GetBinContent(i) + r->GetBinError(i));
                }
                ymax = std::min(std::max(1.15, 1.15 * ymax), 3.0);
                DrawEffFrame(0.0, 1.0, "#DeltaR", 0.0, ymax, eps_single_text);
                DrawUnityLine(0.0, 1.0);
                auto* leg = new TLegend(0.45, 0.68, 0.92, 0.88);
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.03);
                std::vector<std::string> offscale;
                for (int iy = 0; iy < npt; ++iy) {
                    rs[iy]->SetMarkerStyle(ptmark[iy % ptmark.size()]);
                    rs[iy]->SetMarkerColor(ptcol[iy % ptcol.size()]);
                    rs[iy]->SetLineColor(ptcol[iy % ptcol.size()]);
                    rs[iy]->SetLineWidth(2);
                    rs[iy]->Draw("E1 same");
                    leg->AddEntry(rs[iy], pt_label(iy + 1).c_str(), "lp");
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
                leg->Draw();
                DrawHeadline(headline);
                if (!offscale.empty()) {
                    TLatex note; note.SetNDC(); note.SetTextFont(42); note.SetTextSize(0.026);
                    note.SetTextColor(kGray + 3);
                    double y = 0.63;
                    for (size_t i = 0; i < offscale.size(); i += 2) {
                        std::string txt = (i == 0) ? "above scale (arrows): " : "  ";
                        txt += offscale[i];
                        if (i + 1 < offscale.size()) txt += ", " + offscale[i + 1];
                        note.DrawLatex(0.45, y, txt.c_str()); y -= 0.035;
                    }
                }
                SaveCanvas(c, dir4 + "step4_eps_dr_single_zoom_pair_pt_slices.png");
            }

            // pair-eta panels (each subplot = eta bin, each line = pair-pT bin), zoom + full
            struct Rng { TH3D* n; TH3D* d; TH3D* a; TH3D* b; TH3D* p; TH3D* q; std::string tag; double xhi; };
            const std::vector<Rng> rngs = {{h3zn, h3zd, h3zA, h3zB, h3zP, h3zQ, "zoom", 1.0},
                                          {h3fn, h3fd, h3fA, h3fB, h3fP, h3fQ, "full", 5.75}};
            const int ncol = (int)std::ceil(std::sqrt((double)neta));
            const int nrow = (int)std::ceil((double)neta / ncol);
            for (const auto& R : rngs) {
                TCanvas c(("c_s4_pteta_" + R.tag).c_str(), "", 500 * ncol, 450 * nrow);
                c.Divide(ncol, nrow);
                for (int iz = 1; iz <= neta; ++iz) {
                    c.cd(iz);
                    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
                    std::vector<TH1D*> rs; double ymax = 0.;
                    for (int iy = 1; iy <= npt; ++iy) {
                        TH1D* r = cell_ratio(R.n, R.d, R.a, R.b, R.p, R.q, iy, iz, neta,
                                             Form("s4pe_%s_%s_%d_%d", sample.c_str(), R.tag.c_str(), iy, iz));
                        rs.push_back(r);
                        for (int i = 1; i <= r->GetNbinsX(); ++i)
                            ymax = std::max(ymax, r->GetBinContent(i) + r->GetBinError(i));
                    }
                    ymax = std::min(std::max(1.15, 1.15 * ymax), 3.0);
                    DrawEffFrame(0.0, R.xhi, "#DeltaR", 0.0, ymax, eps_single_text);
                    DrawUnityLine(0.0, R.xhi);
                    for (int iy = 0; iy < npt; ++iy) {
                        TH1D* r = rs[iy];
                        r->SetMarkerStyle(ptmark[iy % ptmark.size()]);
                        r->SetMarkerColor(ptcol[iy % ptcol.size()]);
                        r->SetLineColor(ptcol[iy % ptcol.size()]);
                        r->SetLineWidth(2);
                        r->Draw("E1 same");
                        for (int i = 1; i <= r->GetNbinsX(); ++i) {
                            if (r->GetBinContent(i) <= ymax) continue;
                            auto* ar = new TArrow(r->GetBinCenter(i), ymax * 0.88,
                                                  r->GetBinCenter(i), ymax * 0.985, 0.008, "|>");
                            ar->SetLineColor(ptcol[iy % ptcol.size()]);
                            ar->SetFillColor(ptcol[iy % ptcol.size()]); ar->Draw();
                        }
                    }
                    TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.050);
                    tl.DrawLatex(0.17, 0.86, eta_label(iz).c_str());
                    if (iz == 1) {
                        auto* leg = new TLegend(0.42, 0.66, 0.92, 0.90);
                        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
                        for (int iy = 0; iy < npt; ++iy) leg->AddEntry(rs[iy], pt_label(iy + 1).c_str(), "lp");
                        leg->Draw();
                    }
                }
                c.cd(0);
                TLatex st; st.SetNDC(); st.SetTextFont(42); st.SetTextSize(0.016);
                st.DrawLatex(0.03, 0.978, (headline + "  --  " + R.tag + " #DeltaR, " + eps_single_text +
                             " in (p_{T}^{pair}, #eta^{pair}) cells").c_str());
                SaveCanvas(c, dir4 + "step4_eps_dr_single_" + R.tag + "_pair_eta_pt.png");
            }

            // ---- plateau + fluctuation tables (full-dR 3D, window [kPlateauLo,kPlateauHi]) ----
            struct Plat { double mean, err, rms; int nb; };
            auto cell_plateau = [&](int iy, int iz) -> Plat {
                TH1D* r = cell_ratio(h3fn, h3fd, h3fA, h3fB, h3fP, h3fQ, iy, iz, neta,
                                     Form("s4plat_%s_%d_%d", sample.c_str(), iy, iz));
                double sw = 0, swv = 0; std::vector<std::pair<double,double>> vw;
                for (int i = 1; i <= r->GetNbinsX(); ++i) {
                    const double xc = r->GetBinCenter(i);
                    if (xc < kPlateauLo || xc > kPlateauHi) continue;
                    const double v = r->GetBinContent(i), e = r->GetBinError(i);
                    if (e <= 0. || v == 0.) continue;
                    const double w = 1. / (e * e); sw += w; swv += w * v; vw.push_back({v, w});
                }
                delete r;
                if (sw <= 0.) return Plat{ -1, -1, -1, 0 };
                const double mean = swv / sw, err = std::sqrt(1. / sw);
                double swd = 0; for (auto& p : vw) swd += p.second * (p.first - mean) * (p.first - mean);
                return Plat{ mean, err, std::sqrt(swd / sw), (int)vw.size() };
            };
            std::vector<std::vector<Plat>> P(npt, std::vector<Plat>(neta));
            for (int iy = 1; iy <= npt; ++iy)
                for (int iz = 1; iz <= neta; ++iz) P[iy - 1][iz - 1] = cell_plateau(iy, iz);
            auto eta_hdr = [&](int iz){ return Form("[%.1f,%.1f)",
                h3fn->GetZaxis()->GetBinLowEdge(iz), h3fn->GetZaxis()->GetBinUpEdge(iz)); };
            auto pt_hdr  = [&](int iy){ return Form("pTpair[%.0f,%.0f)",
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
                os << "# each cell = stat_err(mean) / rms_scatter / n_bins.\n";
                os << std::left << std::setw(24) << "pair-eta \\ pair-pT";
                for (int iy = 1; iy <= npt; ++iy) os << std::setw(24) << pt_hdr(iy);
                os << "\n";
                for (int iz = 1; iz <= neta; ++iz) {
                    os << std::left << std::setw(24) << eta_hdr(iz);
                    for (int iy = 1; iy <= npt; ++iy) {
                        const Plat& p = P[iy - 1][iz - 1];
                        os << std::setw(24) << (p.nb > 0 ? Form("%.4f/%.4f/%d", p.err, p.rms, p.nb) : "--");
                    }
                    os << "\n";
                }
            };
            write_table_A(std::cout);
            { std::ofstream ofa(dir4 + "step4_plateau_pair_eta_pt.txt");             write_table_A(ofa); }
            { std::ofstream ofb(dir4 + "step4_plateau_fluctuation_pair_eta_pt.txt"); write_table_B(ofb); }
            std::cout << "  wrote " << dir4 << "step4_plateau_pair_eta_pt.txt + "
                      << "step4_plateau_fluctuation_pair_eta_pt.txt\n";
            fmc4->Close();
        }
    }

    fmc->Close(); fmc3->Close(); ffit->Close(); fdata->Close();
    std::cout << "\nplot_mc_trig_eff(" << sample << ") done.\n";
}
