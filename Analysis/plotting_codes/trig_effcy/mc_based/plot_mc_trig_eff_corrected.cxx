// plot_mc_trig_eff_corrected.cxx
// -----------------------------------------------------------------------------------------
// CORRECTED-MC study of the MC trigger efficiency
// (docs/tracking/mc_trigger_efficiency.md, round-7 Autonomy Contract item 5).
//
// Data and MC disagree on the single-muon mu4 efficiency (§3.1, R3: MC above data by ~10-15% in
// the turn-on, more in the forward endcap). There are therefore two ways to build the per-pair
// trigger weight:
//   (1) NOMINAL: data-derived single-muon eps^nc(pT, q.eta)  x  MC-derived (dR, pair pT, pair
//       eta) corrections eps_dR^cross (Step 3) and eps_dR^single (Step 4);
//   (2) CORRECTED MC: weight every MC muon by SF(pT, q.eta) = eps_data/eps_MC so that the MC
//       reproduces the DATA single-muon efficiency, and take everything from that MC.
//
// This macro answers the two questions the choice hinges on.
//
//   Q1 -- do the two give the SAME single-muon efficiency?
//         Step 1 replotted with the CORRECTED MC (numerator weighted by SF, denominator
//         untouched => eps_corr = eps_MC*<SF> ~ eps_data) overlaid on the data tag-and-probe,
//         with the corrected MC's OWN turn-on fit drawn. Quantified in
//         step1_corrected_vs_data.txt: the residual can only come from fit quality, binning and
//         the floor/cap guards, and the file gives all of it in numbers.
//
//   Q2 -- does correcting the MC change the dR CORRECTION TERMS? This is the crucial one: the
//         dR correction is measured in the MC world and applied in the data world, so it is only
//         legitimate if it does not depend on the single-muon normalisation. Algebraically it
//         does not: with the corrected numerator weight w*SF/eps_corr, if eps_corr were exactly
//         eps_data then SF/eps_corr = 1/eps_MC, the ORIGINAL weight, and the corrections would be
//         identical bin by bin. Steps 3 and 4 are replotted as 1 PNG per pair-pT bin, 1 subplot
//         per pair-eta bin, original vs corrected overlaid, with a PULL pad
//         (corrected - original)/sigma_original underneath and a per-cell pull table.
//
// ERROR BARS (round 7): the Step-3/4 ratios use the conditional (binomial-correct) form
// SetConditionalRatioErrors() -- the numerator is a re-weighted SUBSET of the denominator, so
// TH1::Divide's independent propagation over-states them by 1.5-3.2x. The helper is COPIED
// verbatim from plot_mc_trig_eff.cxx (that macro is owned elsewhere; this one must not edit it).
// In corrected mode the errA/errB/covP/covQ histograms are the generalised terms booked by
// FillMCTrigEffHists.cxx (A = sum a^2, B = sum a^2*eps_MC, covP = sum 2a1a2,
// covQ = sum 2a1a2*eps_MC1*eps_MC2), which reduce exactly to the nominal ones at SF = 1.
//
// WP config (registry: Analysis/docs/muon_wp_registry.md): use_tight_wp default TRUE (Tight
// nominal, unsuffixed inputs, plots at the top of each step dir); false selects the _medium_wp
// MC inputs AND the WP-matched _medium_wp DATA reference, into a medium/ subdirectory.
//
// Outputs (NEW directories; the nominal step1/2/3/4 plot dirs are never touched):
//   <out_base>/step1_corrected_mc/[medium/]
//   <out_base>/step3_corrected_mc/[medium/]
//   <out_base>/step4_corrected_mc/[medium/]
//
// Compile/run (ACLiC, from this directory):
//   root -l -b -q 'plot_mc_trig_eff_corrected.cxx+("pp_full")'
//   root -l -b -q 'plot_mc_trig_eff_corrected.cxx+("overlay", false)'
// -----------------------------------------------------------------------------------------

#include <TArrow.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// ================================================================= generic helpers
// (BayesEff / StyleGraph / DrawEffFrame / SplitPadForRatio / DrawRatioFrame /
//  SetConditionalRatioErrors / DivideGraphClean / MarkOffScale / DrawHeadline /
//  PlateauWeightedMean are the plot_mc_trig_eff.cxx conventions, reproduced here so that this
//  macro is self-contained and that macro is left untouched.)

template <typename T>
T* GetObj(TFile* f, const std::string& name)
{
    T* obj = dynamic_cast<T*>(f->Get(name.c_str()));
    if (!obj)
        throw std::runtime_error("plot_mc_trig_eff_corrected: missing object '" + name +
                                 "' in file " + f->GetName());
    return obj;
}

TFile* OpenFile(const std::string& path)
{
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie())
        throw std::runtime_error("plot_mc_trig_eff_corrected: cannot open " + path);
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

std::pair<TPad*, TPad*> SplitPadForRatio(const std::string& tag, bool logx, double split = 0.32)
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

// Frame for a ratio/pull pad: labels/titles scaled by 1/split so they match the main pad.
// `ref` is the horizontal reference line (1 for a ratio, 0 for a pull).
TH1* DrawRatioFrame(double xlo, double xhi, const std::string& xtitle,
                    const std::string& ytitle, double ylo, double yhi,
                    double ref = 1.0, double split = 0.32)
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
    auto* l = new TLine(xlo, ref, xhi, ref);
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->Draw("same");
    return fr;
}

// =============================================================================
// CONDITIONAL (binomial-correct) ERROR ON AN INVERSE-WEIGHTED EFFICIENCY RATIO  (round 7)
// Verbatim from plot_mc_trig_eff.cxx (see that file's header for the full derivation):
//     Var = A - R*B (+ covP - R^2*covQ for Step 4),  e_R = sqrt(Var)/D
// with the k=n boundary fallback to the "1/n rule" on the effective denominator count.
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

        const double scale = a + R * b + std::fabs(cp) + R * R * cq;
        if (var > 1e-6 * scale) { r->SetBinError(i, std::sqrt(var) / D); continue; }

        const double eD = den->GetBinError(i);
        const double neff = (eD > 0.) ? (D / eD) * (D / eD) : 1.0;
        r->SetBinError(i, (neff > 0.) ? R / neff : 0.);
    }
}

// Point-by-point ratio of two efficiency graphs, MATCHED BY X (TGraphAsymmErrors::Divide skips
// empty-denominator bins, so indices do NOT correspond). Verbatim from plot_mc_trig_eff.cxx.
TGraphAsymmErrors* DivideGraphClean(TGraphAsymmErrors* gn, TGraphAsymmErrors* gd)
{
    auto* gr = new TGraphAsymmErrors();
    int k = 0;
    for (int i = 0; i < gn->GetN(); ++i) {
        double xn, yn;
        gn->GetPoint(i, xn, yn);
        int j = -1;
        for (int m = 0; m < gd->GetN(); ++m) {
            double xd, yd;
            gd->GetPoint(m, xd, yd);
            if (std::fabs(xn - xd) <= 1e-6 * std::max(1.0, std::fabs(xd))) { j = m; break; }
        }
        if (j < 0) continue;
        double xd, yd;
        gd->GetPoint(j, xd, yd);
        if (yd <= 0.) continue;
        const double en = 0.5 * (gn->GetErrorYhigh(i) + gn->GetErrorYlow(i));
        const double ed = 0.5 * (gd->GetErrorYhigh(j) + gd->GetErrorYlow(j));
        const double r  = yn / yd;
        const double rel_n = (yn > 0.) ? (en / yn) : 0.;
        const double er = (yn > 0.) ? r * std::sqrt(rel_n * rel_n + (ed / yd) * (ed / yd))
                                    : en / yd;
        gr->SetPoint(k, xn, r);
        gr->SetPointError(k, gn->GetErrorXlow(i), gn->GetErrorXhigh(i), er, er);
        ++k;
    }
    return gr;
}

// A zoomed frame hides points outside it. Mark every such point with an arrow at the frame edge:
// silently dropping data from a physics figure is not allowed.
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

// Same, for a histogram; returns the "x = value +- error" strings of the marked points so the
// caller can list them on the canvas.
std::vector<std::string> MarkOffScaleHist(TH1D* h, double ylo, double yhi, Color_t col,
                                          const char* fmt = "#DeltaR=%.2f: %.2f #pm %.2f")
{
    std::vector<std::string> out;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        const double v = h->GetBinContent(i);
        if (h->GetBinError(i) == 0. && v == 0.) continue;   // empty bin: nothing measured
        if (v <= yhi && v >= ylo) continue;
        const double x = h->GetBinCenter(i);
        const bool up = (v > yhi);
        const double y0 = up ? ylo + 0.88 * (yhi - ylo) : ylo + 0.12 * (yhi - ylo);
        const double y1 = up ? ylo + 0.98 * (yhi - ylo) : ylo + 0.02 * (yhi - ylo);
        auto* ar = new TArrow(x, y0, x, y1, 0.008, "|>");
        ar->SetLineColor(col);
        ar->SetFillColor(col);
        ar->SetLineWidth(2);
        ar->Draw();
        out.push_back(Form(fmt, x, v, h->GetBinError(i)));
    }
    return out;
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
    c.SaveAs(path.c_str());     // PNG only (project convention: never also PDF)
    std::cout << "  wrote " << path << std::endl;
}

size_t GlyphLength(const std::string& s)
{
    size_t n = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '{' || s[i] == '}' || s[i] == '^' || s[i] == '_') continue;
        if (s[i] == '#') {
            ++i;
            while (i < s.size() && std::isalpha(static_cast<unsigned char>(s[i]))) ++i;
            --i;
        }
        ++n;
    }
    return n;
}

void DrawHeadline(const std::string& text, double x = 0.12, double y = 0.955, double size = 0.038)
{
    constexpr size_t kFitsAt = 34;
    const size_t len = GlyphLength(text);
    if (len > kFitsAt) size *= static_cast<double>(kFitsAt) / static_cast<double>(len);
    size = std::max(size, 0.030);
    TLatex tl;
    tl.SetNDC();
    tl.SetTextSize(size);
    tl.SetTextFont(42);
    tl.DrawLatex(x, y, text.c_str());
}

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
    if (sumw <= 0.) return {-1., -1.};
    return {sumwv / sumw, std::sqrt(1. / sumw)};
}

// ================================================================= config

struct SampleCfg {
    std::string mc_dir;
    std::string mc_label;
    std::string data_hist_file;   // T&P hists + graphs (num/denom, g_..._divided)
    std::string data_fit_file;    // T&P TF1s (f_..._divided)
    std::string ctr;              // "" or "_ctr0_5"
    std::string out_base;
    std::string sample_text;
    std::string data_text;
    std::string eps_dr_text;      // eps_dR symbol of Step 3 for this sample
};

// The DATA reference must be at the SAME working point as the MC (§3.0(d)) -- both data files
// are WP-keyed. Paths are the analysis's own (RDFBasedHistFillingPP.cxx:306/323,
// RDFBasedHistFillingPbPb.cxx:678/694) and match plot_mc_trig_eff.cxx MakeCfg.
SampleCfg MakeCfg(const std::string& sample, bool use_tight_wp)
{
    const std::string wp = use_tight_wp ? "" : "_medium_wp";
    SampleCfg c;
    if (sample == "pp_full" || sample == "pp") {
        const bool full = (sample == "pp_full");
        c.mc_dir   = full ? "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/"
                          : "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/";
        c.mc_label = full ? "pp24_full" : "pp24";
        c.data_hist_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
                           "histograms_real_pairs_pp_2024_single_mu4_fine_q_eta_bin" + wp + ".root";
        c.data_fit_file  = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
                           "trg_effcy_pT_fitting_to_erf_plus_log/single_mu_effcy_pT_fit" + wp + ".root";
        c.ctr         = "";
        c.out_base    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                        "pp_trigger_efficiency/mc_based/";
        c.sample_text = full ? "Pythia8 pp24 fullsim (FULL sample)" : "Pythia8 pp24 fullsim";
        c.data_text   = "pp24 data";
        c.eps_dr_text = "#varepsilon_{#DeltaR}^{2mu4}";
    } else if (sample == "overlay") {
        c.mc_dir   = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/";
        c.mc_label = "hijing_overlay_pbpb23";
        c.data_hist_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/"
                           "histograms_real_pairs_pbpb_2023_single_mu4_fine_q_eta_bin" + wp + ".root";
        c.data_fit_file  = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/"
                           "trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit" + wp + ".root";
        c.ctr         = "_ctr0_5";   // D2: the overlay compares ONLY to PbPb23 data 0-5%
        c.out_base    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                        "pbpb_trigger_efficiency/mc_based/";
        c.sample_text = "HIJING overlay Pb+Pb23 cond., 0-5%";
        c.data_text   = "Pb+Pb23 data 0-5%";
        c.eps_dr_text = "#varepsilon_{#DeltaR}^{cross}";
    } else if (sample == "noovl") {
        c.mc_dir   = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_no_overlay_test_sample/";
        c.mc_label = "r17663_no_overlay";
        c.data_hist_file = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
                           "histograms_real_pairs_pp_2024_single_mu4_fine_q_eta_bin" + wp + ".root";
        c.data_fit_file  = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
                           "trg_effcy_pT_fitting_to_erf_plus_log/single_mu_effcy_pT_fit" + wp + ".root";
        c.ctr         = "";
        c.out_base    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/"
                        "r17663_no_overlay_trigger_efficiency/mc_based/";
        c.sample_text = "Pythia8 pp, r17663 (HI cond., no overlay)";
        c.data_text   = "pp24 data";
        c.eps_dr_text = "#varepsilon_{#DeltaR}^{2mu4}";
    } else {
        throw std::runtime_error("plot_mc_trig_eff_corrected: sample must be 'pp', 'pp_full', "
                                 "'overlay' or 'noovl', got " + sample);
    }
    return c;
}

// charge <-> data-sign mapping: sign1 = mu+, sign2 = mu- (data sign convention)
const std::vector<std::string> kCharges   = {"muplus", "muminus"};
const std::vector<std::string> kDataSigns = {"sign1", "sign2"};
const std::vector<std::string> kChargeTex = {"#mu^{+}", "#mu^{-}"};

// fine q.eta bins (CommonEffcyConfig.h q_eta_proj_ranges_fine_excl_gap); the forward
// (-2.4,-2.0) bin is split into two (round-5 change #2) -> 11 bins.
const std::vector<std::string> kQEtaSuffix = {
    "minus2_40_TO_minus2_20", "minus2_20_TO_minus2_00",
    "minus2_00_TO_minus1_60", "minus1_60_TO_minus1_30",
    "minus0_90_TO_minus0_50", "minus0_50_TO_minus0_10", "0_10_TO_0_50",
    "0_50_TO_1_00", "1_30_TO_1_60", "1_60_TO_2_00", "2_00_TO_2_20"};
const std::vector<std::pair<double,double>> kQEtaRange = {
    {-2.4,-2.2},{-2.2,-2.0},{-2.0,-1.6},{-1.6,-1.3},{-0.9,-0.5},{-0.5,-0.1},
    { 0.1, 0.5},{ 0.5, 1.0},{ 1.3, 1.6},{ 1.6, 2.0},{ 2.0, 2.2}};

const Color_t kCorrColor = kRed + 1;     // corrected MC
const Color_t kOrigColor = kBlue + 1;    // original (uncorrected) MC
const Color_t kDataColor = kBlack;       // data tag-and-probe

// Step-3/4 plateau window (§3.3 diagnostic 1), same as the nominal macro.
const double kPlateauLo = 1.0;
const double kPlateauHi = 4.0;

// ---- Step-3/4 cell projection: eps(dR) in one (pair-pT bin iy, pair-eta bin iz) cell --------
// Central value = num/denom; ERROR = conditional/binomial-correct form. covP/covQ are Step-4
// only (both legs of a pair share a dR bin); pass nullptr for Step 3.
TH1D* CellRatio(TH3D* hn, TH3D* hd, TH3D* ha, TH3D* hb, TH3D* hp, TH3D* hq,
                int iy, int iz, const std::string& nm)
{
    TH1D* n = hn->ProjectionX((nm + "_n").c_str(), iy, iy, iz, iz, "e");
    TH1D* d = hd->ProjectionX((nm + "_d").c_str(), iy, iy, iz, iz, "e");
    TH1D* a = ha->ProjectionX((nm + "_a").c_str(), iy, iy, iz, iz, "e");
    TH1D* b = hb->ProjectionX((nm + "_b").c_str(), iy, iy, iz, iz, "e");
    TH1D* p = hp ? hp->ProjectionX((nm + "_p").c_str(), iy, iy, iz, iz, "e") : nullptr;
    TH1D* q = hq ? hq->ProjectionX((nm + "_q").c_str(), iy, iy, iz, iz, "e") : nullptr;
    auto* r = (TH1D*)n->Clone(nm.c_str());
    r->SetDirectory(nullptr);
    r->Divide(d);
    SetConditionalRatioErrors(r, d, a, b, p, q);
    delete n; delete d; delete a; delete b;
    if (p) delete p;
    if (q) delete q;
    return r;
}

// A bundle of the 3D inputs of one step / one dR range.
struct Step34Hists {
    TH3D *num, *den, *errA, *errB, *covP, *covQ;
};
Step34Hists LoadStep34(TFile* f, const std::string& base, const std::string& range, bool with_cov)
{
    const std::string p = base + range + "_vs_pt_eta_";
    Step34Hists h;
    h.num  = GetObj<TH3D>(f, p + "num");
    h.den  = GetObj<TH3D>(f, p + "denom");
    h.errA = GetObj<TH3D>(f, p + "errA");
    h.errB = GetObj<TH3D>(f, p + "errB");
    h.covP = with_cov ? GetObj<TH3D>(f, p + "covP") : nullptr;
    h.covQ = with_cov ? GetObj<TH3D>(f, p + "covQ") : nullptr;
    return h;
}

} // namespace

// ================================================================= main

void plot_mc_trig_eff_corrected(const std::string& sample = "pp_full", bool use_tight_wp = true)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gErrorIgnoreLevel = kWarning;

    const SampleCfg cfg = MakeCfg(sample, use_tight_wp);
    const std::string wp_suf  = use_tight_wp ? "" : "_medium_wp";
    const std::string wp_text = use_tight_wp ? "Tight muons" : "Medium muons";
    const std::string wp_dir  = use_tight_wp ? "" : "medium/";
    const std::string headline = cfg.sample_text + ", " + wp_text;

    // ---- inputs ----
    TFile* fcorr1 = OpenFile(cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf + "_corrected.root");
    TFile* fnom1  = OpenFile(cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf + ".root");
    TFile* fcfit  = OpenFile(cfg.mc_dir + "single_mu_effcy_pT_fit_mc_corrected" + wp_suf + ".root");
    TFile* fnfit  = OpenFile(cfg.mc_dir + "single_mu_effcy_pT_fit_mc" + wp_suf + ".root");
    TFile* fdata  = OpenFile(cfg.data_hist_file);
    TFile* fdfit  = OpenFile(cfg.data_fit_file);

    const std::string dir1 = cfg.out_base + "step1_corrected_mc/" + wp_dir;
    const std::string dir3 = cfg.out_base + "step3_corrected_mc/" + wp_dir;
    const std::string dir4 = cfg.out_base + "step4_corrected_mc/" + wp_dir;
    for (const auto& d : {dir1, dir3, dir4}) gSystem->mkdir(d.c_str(), kTRUE);

    const std::string corr_leg = "corrected MC (#times #varepsilon_{data}/#varepsilon_{MC})";
    const std::string orig_leg = "original MC, P(mu4 | reco #mu)";
    const std::string data_leg = cfg.data_text + " T&P P(2mu4 | mu4 tag, #DeltaR>0.8)";

    // ==================================================================================
    // Q1 / Step 1: does the corrected MC reproduce the DATA single-muon efficiency?
    // ==================================================================================
    std::cout << "\n===== Q1: Step 1 corrected MC vs data (" << sample << ", " << wp_text
              << ") =====\n";

    struct Var { std::string mc, data, xtitle; double xlo, xhi; bool logx; };
    const std::vector<Var> vars = {
        {"pt",  "pt2nd",  "p_{T} [GeV]", 4.0, 60.0, true},
        {"eta", "eta2nd", "#eta",       -2.4,  2.4, false},
        {"phi", "phi2nd", "#phi",       -M_PI, M_PI, false}};

    // integrated corrected/data and original/data, per charge, for the summary table
    std::vector<std::array<double,4>> integ(2);   // {eps_corr, eps_orig, eps_data, n/a}

    for (const auto& v : vars) {
        TCanvas c(("c_s1c_" + v.mc).c_str(), "", 1400, 800);
        c.Divide(2, 1);
        for (int ic = 0; ic < 2; ++ic) {
            c.cd(ic + 1);
            auto pads = SplitPadForRatio("s1c_" + v.mc + "_" + std::to_string(ic), v.logx);

            TH1D* cnum = GetObj<TH1D>(fcorr1, "h_mc_" + v.mc + "_num_"   + kCharges[ic]);
            TH1D* cden = GetObj<TH1D>(fcorr1, "h_mc_" + v.mc + "_denom_" + kCharges[ic]);
            TH1D* onum = GetObj<TH1D>(fnom1,  "h_mc_" + v.mc + "_num_"   + kCharges[ic]);
            TH1D* oden = GetObj<TH1D>(fnom1,  "h_mc_" + v.mc + "_denom_" + kCharges[ic]);
            TH1D* dnum = GetObj<TH1D>(fdata, "h_" + v.data + cfg.ctr + "_" + kDataSigns[ic] + "_2mu4_sepr");
            TH1D* dden = GetObj<TH1D>(fdata, "h_" + v.data + cfg.ctr + "_" + kDataSigns[ic] + "_mu4_sepr");

            // The CORRECTED numerator is re-weighted by SF and can EXCEED the denominator in a
            // bin where SF > 1, so a Bayesian/binomial divide is undefined (BayesDivide rejects
            // the pair and returns an empty graph). Build the estimator under study explicitly,
            // eps_corr = sum_fired w*SF / sum_all w, with the CONDITIONAL error
            // e = sqrt(A - B)/D, A = sum_fired (w*SF)^2, B = sum_fired (w*SF)^2 * eps_MC
            // (booked by FillMCTrigEffHists in corrected mode; same form as the fitter's
            // CorrectedEffGraph, incl. the "1/n rule" fallback for BOTH binomial boundaries --
            // max(eff, 1)/n_eff, so that an eff = 0 bin gets a finite error instead of 0).
            TH1D* cA = GetObj<TH1D>(fcorr1, "h_mc_" + v.mc + "_errA_" + kCharges[ic]);
            TH1D* cB = GetObj<TH1D>(fcorr1, "h_mc_" + v.mc + "_errB_" + kCharges[ic]);
            auto* gcorr = new TGraphAsymmErrors();
            for (int i = 1, k = 0; i <= cden->GetNbinsX(); ++i) {
                const double D = cden->GetBinContent(i);
                const double bw = cden->GetBinWidth(i);
                if (D <= 0. || bw <= 0.) continue;
                const double N = cnum->GetBinContent(i);
                const double a = cA->GetBinContent(i), b = cB->GetBinContent(i);
                const double var = a - b;
                double e;
                if (var > 1e-6 * (a + b)) {
                    e = std::sqrt(var) / D;
                } else {
                    const double eD = cden->GetBinError(i);
                    const double neff = (eD > 0.) ? (D / eD) * (D / eD) : 1.0;
                    e = (neff > 0.) ? std::max(N / D, 1.0) / neff : 0.;
                }
                gcorr->SetPoint(k, cden->GetBinCenter(i), N / D);
                gcorr->SetPointError(k, 0.5 * bw, 0.5 * bw, e, e);
                ++k;
            }
            auto* gorig = BayesEff(onum, oden);
            auto* gdata = BayesEff(dnum, dden);
            StyleGraph(gcorr, kCorrColor, 21);
            StyleGraph(gorig, kOrigColor, 22, 0.8);
            StyleGraph(gdata, kDataColor, 20);

            if (v.mc == "pt") {
                integ[ic] = {cden->Integral() > 0 ? cnum->Integral() / cden->Integral() : -1.,
                             oden->Integral() > 0 ? onum->Integral() / oden->Integral() : -1.,
                             dden->Integral() > 0 ? dnum->Integral() / dden->Integral() : -1.,
                             0.};
            }

            pads.first->cd();
            DrawEffFrame(v.xlo, v.xhi, "", 0.0, 1.1, "efficiency", true);
            DrawUnityLine(v.xlo, v.xhi);
            gdata->Draw("PZ same");
            gorig->Draw("PZ same");
            gcorr->Draw("PZ same");

            DrawHeadline(headline + ", " + kChargeTex[ic], 0.14, 0.955, 0.05);
            auto* leg = new TLegend(0.18, 0.08, 0.93, 0.32);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.034);
            leg->AddEntry(gdata, data_leg.c_str(), "lp");
            leg->AddEntry(gcorr, corr_leg.c_str(), "lp");
            leg->AddEntry(gorig, orig_leg.c_str(), "lp");
            leg->Draw();

            pads.second->cd();
            DrawRatioFrame(v.xlo, v.xhi, v.xtitle, "MC / data", 0.5, 2.2);
            auto* grc = DivideGraphClean(gcorr, gdata);
            auto* gro = DivideGraphClean(gorig, gdata);
            StyleGraph(gro, kOrigColor, 22, 0.8);
            StyleGraph(grc, kCorrColor, 21, 0.9);
            gro->Draw("PZ same");
            grc->Draw("PZ same");
            MarkOffScale(gro, 0.5, 2.2, kOrigColor);
            MarkOffScale(grc, 0.5, 2.2, kCorrColor);
            c.cd(ic + 1);
        }
        SaveCanvas(c, dir1 + "step1_corrected_eff_" + v.mc + ".png");
    }

    // --- pT turn-on in the 11 fine q.eta bins, per charge (3x4 grid: 11 panels + legend) ---
    for (int ic = 0; ic < 2; ++ic) {
        TCanvas c(("c_s1c_qeta_" + kCharges[ic]).c_str(), "", 1500, 2100);
        c.Divide(3, 4);   // ncols = 3, nrows = 4 (nrows >= ncols)
        for (size_t iq = 0; iq < kQEtaSuffix.size(); ++iq) {
            c.cd(static_cast<int>(iq) + 1);
            auto pads = SplitPadForRatio("s1cq_" + kCharges[ic] + "_" + std::to_string(iq), true);

            auto* gc = (TGraphAsymmErrors*)GetObj<TGraphAsymmErrors>(
                fcfit, "g_mc_pt_vs_q_eta_" + kCharges[ic] + "_" + kQEtaSuffix[iq])->Clone();
            auto* go = (TGraphAsymmErrors*)GetObj<TGraphAsymmErrors>(
                fnfit, "g_mc_pt_vs_q_eta_" + kCharges[ic] + "_" + kQEtaSuffix[iq])->Clone();
            auto* gd = (TGraphAsymmErrors*)GetObj<TGraphAsymmErrors>(
                fdata, "g_pt2nd_vs_q_eta2nd" + cfg.ctr + "_" + kDataSigns[ic] +
                       "_2mu4_sepr_py_" + kQEtaSuffix[iq] + "_divided")->Clone();
            auto* fc = (TF1*)GetObj<TF1>(fcfit, "f_mc_pt_vs_q_eta_" + kCharges[ic] + "_" +
                                                kQEtaSuffix[iq])->Clone();
            auto* fd = (TF1*)GetObj<TF1>(fdfit, "f_pt2nd_vs_q_eta2nd" + cfg.ctr + "_" +
                                                kDataSigns[ic] + "_2mu4_sepr_py_" +
                                                kQEtaSuffix[iq] + "_divided")->Clone();
            // The graphs carry their own stored fit (they were written after TGraph::Fit) and
            // ROOT would auto-draw it with "PZ same" -- a second, unrequested red curve. Strip it
            // and draw only the curves we mean to show.
            gc->GetListOfFunctions()->Clear();
            go->GetListOfFunctions()->Clear();
            gd->GetListOfFunctions()->Clear();
            StyleGraph(gc, kCorrColor, 21, 0.7);
            StyleGraph(go, kOrigColor, 22, 0.7);
            StyleGraph(gd, kDataColor, 20, 0.7);

            pads.first->cd();
            gPad->SetLogx();
            DrawEffFrame(4.0, 60.0, "", 0.0, 1.1, "efficiency", true);
            DrawUnityLine(4.0, 60.0);
            fd->SetLineColor(kDataColor); fd->SetLineWidth(1); fd->SetLineStyle(2);
            fd->Draw("same");
            fc->SetLineColor(kCorrColor); fc->SetLineWidth(2);
            fc->Draw("same");
            gd->Draw("PZ same");
            go->Draw("PZ same");
            gc->Draw("PZ same");

            TLatex tl;
            tl.SetNDC(); tl.SetTextSize(0.070); tl.SetTextFont(42);
            tl.DrawLatex(0.28, 0.12, Form("%.1f < q#upoint#eta < %.1f",
                                          kQEtaRange[iq].first, kQEtaRange[iq].second));

            pads.second->cd();
            gPad->SetLogx();
            DrawRatioFrame(4.0, 60.0, "p_{T} [GeV]", "MC / data", 0.5, 2.6);
            auto* grc = DivideGraphClean(gc, gd);
            auto* gro = DivideGraphClean(go, gd);
            StyleGraph(gro, kOrigColor, 22, 0.7);
            StyleGraph(grc, kCorrColor, 21, 0.7);
            gro->Draw("PZ same");
            grc->Draw("PZ same");
            MarkOffScale(gro, 0.5, 2.6, kOrigColor);
            MarkOffScale(grc, 0.5, 2.6, kCorrColor);
            c.cd(static_cast<int>(iq) + 1);
        }
        c.cd(12);
        DrawHeadline(headline + ", " + kChargeTex[ic], 0.02, 0.90, 0.046);
        auto* lgc = new TGraphAsymmErrors(); StyleGraph(lgc, kCorrColor, 21);
        auto* lgo = new TGraphAsymmErrors(); StyleGraph(lgo, kOrigColor, 22);
        auto* lgd = new TGraphAsymmErrors(); StyleGraph(lgd, kDataColor, 20);
        auto* leg = new TLegend(0.02, 0.56, 0.98, 0.86);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.042);
        leg->AddEntry(lgd, (cfg.data_text + " T&P (+ fit, dashed)").c_str(), "lp");
        leg->AddEntry(lgc, "corrected MC (+ its own fit, solid)", "lp");
        leg->AddEntry(lgo, "original MC (points only)", "lp");
        leg->Draw();
        TLatex note;
        note.SetNDC(); note.SetTextSize(0.038); note.SetTextFont(42);
        note.DrawLatex(0.02, 0.48, "Corrected MC: every FIRED muon weighted by");
        note.DrawLatex(0.02, 0.42, "SF(p_{T},q#upoint#eta) = #varepsilon_{data}/#varepsilon_{MC}"
                                   " (denominator untouched)");
        note.DrawLatex(0.02, 0.34, "#Rightarrow #varepsilon_{corr} = #varepsilon_{MC}"
                                   "#LTSF#GT #approx #varepsilon_{data}.");
        note.DrawLatex(0.02, 0.26, "Ratio pad: MC / data, both series.");
        note.DrawLatex(0.02, 0.18, "Residual #ne 0 comes from fit quality, binning");
        note.DrawLatex(0.02, 0.12, "and the #varepsilon floor/cap + q#upoint#eta-gap fallbacks.");
        SaveCanvas(c, dir1 + "step1_corrected_eff_pt_in_q_eta_bins_" + kCharges[ic] + ".png");
    }

    // --- Q1 quantification: eps_corr / eps_data at the FIT level (that is what enters Steps 3/4)
    {
        std::ofstream os(dir1 + "step1_corrected_vs_data.txt");
        os << "# Q1: does correcting the MC by SF = eps_data/eps_MC reproduce the DATA "
              "single-muon mu4 efficiency?\n";
        os << "# sample=" << sample << "  WP=" << wp_text << "  " << cfg.sample_text << "\n";
        os << "# (a) INTEGRATED efficiency over the whole selected sample (pT > 4, |eta| < 2.4).\n";
        os << "#     CAUTION -- corr/data is NOT expected to be 1 here, and a deviation is not a\n"
              "#     failure of the correction: eps_corr is the MC-spectrum-weighted average of\n"
              "#     eps_data(pT, q.eta) while eps_data is the DATA-spectrum-weighted average, and\n"
              "#     the two samples have different (pT, q.eta) mixtures. On top of that the data\n"
              "#     denominator cannot be truth-matched and still contains fakes/hadronic muons of\n"
              "#     ~zero trigger efficiency (the documented §3.0 asymmetry), which dilutes\n"
              "#     eps_data. The DIFFERENTIAL comparison in (b) is the meaningful test -- and it\n"
              "#     is also the only one that matters downstream, since Steps 3/4 use eps at the\n"
              "#     muon's own (pT, q.eta), never an integrated value.\n";
        os << std::left << std::setw(10) << "charge" << std::setw(14) << "eps_corr"
           << std::setw(14) << "eps_orig" << std::setw(14) << "eps_data"
           << std::setw(14) << "corr/data" << std::setw(14) << "orig/data" << "\n";
        for (int ic = 0; ic < 2; ++ic)
            os << std::left << std::setw(10) << kCharges[ic]
               << std::setw(14) << Form("%.4f", integ[ic][0])
               << std::setw(14) << Form("%.4f", integ[ic][1])
               << std::setw(14) << Form("%.4f", integ[ic][2])
               << std::setw(14) << Form("%.4f", integ[ic][2] > 0 ? integ[ic][0] / integ[ic][2] : -1.)
               << std::setw(14) << Form("%.4f", integ[ic][2] > 0 ? integ[ic][1] / integ[ic][2] : -1.)
               << "\n";

        os << "\n# (b) FIT-LEVEL ratio eps_corr(pT)/eps_data(pT) per fine q.eta bin, sampled on a\n"
              "#     100-point log grid over pT in [4,60] GeV (the range both fits are defined on).\n"
              "#     This is the quantity that propagates into Steps 3/4: the corrected inverse\n"
              "#     weight is SF/eps_corr = (eps_data/eps_MC)/eps_corr, which equals the ORIGINAL\n"
              "#     1/eps_MC exactly when this ratio is 1.\n";
        os << "#     Columns: mean ratio; max |ratio-1| over the WHOLE range; max |ratio-1|\n"
              "#     restricted to eps_data > 0.2, i.e. excluding the first ~1 GeV of the turn-on\n"
              "#     where eps_data is a few percent and any small ABSOLUTE difference is a huge\n"
              "#     RELATIVE one (that region is also where the two functional forms differ most\n"
              "#     and where the eps floor bites); and the pT of the unrestricted maximum.\n";
        os << std::left << std::setw(26) << "q.eta bin" << std::setw(10) << "charge"
           << std::setw(14) << "<corr/data>" << std::setw(14) << "max|dev|"
           << std::setw(16) << "max|dev|(e>0.2)"
           << std::setw(12) << "at pT" << std::setw(14) << "<orig/data>"
           << std::setw(14) << "max|dev|orig" << "\n";

        double worst_corr = 0., worst_pt = 0.; std::string worst_bin;
        double worst_corr_hi = 0.; std::string worst_bin_hi;
        double sum_corr = 0.; int nsum = 0;
        for (size_t iq = 0; iq < kQEtaSuffix.size(); ++iq) {
            for (int ic = 0; ic < 2; ++ic) {
                auto* fc = GetObj<TF1>(fcfit, "f_mc_pt_vs_q_eta_" + kCharges[ic] + "_" + kQEtaSuffix[iq]);
                auto* fo = GetObj<TF1>(fnfit, "f_mc_pt_vs_q_eta_" + kCharges[ic] + "_" + kQEtaSuffix[iq]);
                auto* fd = GetObj<TF1>(fdfit, "f_pt2nd_vs_q_eta2nd" + cfg.ctr + "_" + kDataSigns[ic] +
                                              "_2mu4_sepr_py_" + kQEtaSuffix[iq] + "_divided");
                double sc = 0., so = 0., mc_ = 0., mchi = 0., mo = 0., ptmax = 0.;
                const int NP = 100;
                for (int i = 0; i < NP; ++i) {
                    // log grid, and each TF1 evaluated CONTINUOUSLY at the exact pT
                    const double pt = 4.0 * std::pow(60.0 / 4.0, double(i) / (NP - 1));
                    const double ed = std::min(1.0, fd->Eval(pt));
                    if (ed <= 0.) continue;
                    const double rc = std::min(1.0, fc->Eval(pt)) / ed;
                    const double ro = std::min(1.0, fo->Eval(pt)) / ed;
                    sc += rc; so += ro;
                    if (std::fabs(rc - 1.) > mc_) { mc_ = std::fabs(rc - 1.); ptmax = pt; }
                    if (ed > 0.2) mchi = std::max(mchi, std::fabs(rc - 1.));
                    mo = std::max(mo, std::fabs(ro - 1.));
                }
                sc /= NP; so /= NP;
                sum_corr += std::fabs(sc - 1.); ++nsum;
                if (mc_ > worst_corr) { worst_corr = mc_; worst_pt = ptmax;
                                        worst_bin = kQEtaSuffix[iq] + " " + kCharges[ic]; }
                if (mchi > worst_corr_hi) { worst_corr_hi = mchi;
                                            worst_bin_hi = kQEtaSuffix[iq] + " " + kCharges[ic]; }
                os << std::left << std::setw(26) << Form("[%.1f,%.1f)", kQEtaRange[iq].first,
                                                         kQEtaRange[iq].second)
                   << std::setw(10) << kCharges[ic]
                   << std::setw(14) << Form("%.4f", sc)
                   << std::setw(14) << Form("%.4f", mc_)
                   << std::setw(16) << Form("%.4f", mchi)
                   << std::setw(12) << Form("%.1f", ptmax)
                   << std::setw(14) << Form("%.4f", so)
                   << std::setw(14) << Form("%.4f", mo) << "\n";
            }
        }
        os << "\n# VERDICT: mean |<corr/data> - 1| over the " << nsum
           << " (q.eta, charge) fits = " << Form("%.4f", nsum ? sum_corr / nsum : -1.)
           << " ; worst point-wise |corr/data - 1| = " << Form("%.4f", worst_corr)
           << " at pT = " << Form("%.1f", worst_pt) << " GeV in " << worst_bin
           << " ; worst point-wise where eps_data > 0.2 = " << Form("%.4f", worst_corr_hi)
           << " in " << worst_bin_hi << "\n";
        std::cout << "  Q1 fit-level: mean |<eps_corr/eps_data> - 1| = "
                  << (nsum ? sum_corr / nsum : -1.) << " , worst point-wise = " << worst_corr
                  << " (" << worst_bin << ", pT " << worst_pt << " GeV)"
                  << " , worst where eps_data>0.2 = " << worst_corr_hi
                  << " (" << worst_bin_hi << ")\n";
        std::cout << "  wrote " << dir1 << "step1_corrected_vs_data.txt\n";
    }

    // ==================================================================================
    // Q2 / Steps 3 and 4: does correcting the MC change the dR CORRECTION TERMS?
    // 1 PNG per pair-pT bin, 1 subplot per pair-eta bin, original vs corrected overlaid,
    // with a PULL pad (corrected - original)/sigma_original underneath.
    // ==================================================================================
    struct StepDef {
        std::string key;        // "step3" | "step4"
        std::string base;       // hist base name
        std::string file_tag;   // "_step3" | "_step4"
        std::string dir;
        std::string ytitle;
        bool with_cov;          // Step 4 has the leg-leg covariance terms
    };
    const std::vector<StepDef> steps = {
        {"step3", "h_mc_dr_",        "_step3", dir3, cfg.eps_dr_text,                 false},
        {"step4", "h_mc_single_dr_", "_step4", dir4, "#varepsilon_{#DeltaR}^{single}", true}};

    for (const auto& S : steps) {
        std::cout << "\n===== Q2: " << S.key << " corrected vs original (" << sample << ", "
                  << wp_text << ") =====\n";
        TFile* fo = OpenFile(cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf +
                             S.file_tag + ".root");
        TFile* fc = OpenFile(cfg.mc_dir + "mc_trig_eff_hists_" + cfg.mc_label + wp_suf +
                             "_corrected" + S.file_tag + ".root");

        struct Rng { std::string tag; double xhi; };
        const std::vector<Rng> rngs = {{"zoom", 1.0}, {"full", 5.75}};

        // pull bookkeeping over ALL cells and BOTH dR ranges, for the summary table
        struct CellStat { double maxpull, medpull, maxabs, plat_o, plat_c; int nb; };
        std::vector<std::vector<CellStat>> stat_full;   // [ipt][ieta], full-dR range

        for (const auto& R : rngs) {
            Step34Hists HO = LoadStep34(fo, S.base, R.tag, S.with_cov);
            Step34Hists HC = LoadStep34(fc, S.base, R.tag, S.with_cov);
            const int npt  = HO.num->GetYaxis()->GetNbins();
            const int neta = HO.num->GetZaxis()->GetNbins();
            auto pt_label  = [&](int iy){ return std::string(Form("%.0f < p_{T}^{pair} < %.0f GeV",
                HO.num->GetYaxis()->GetBinLowEdge(iy), HO.num->GetYaxis()->GetBinUpEdge(iy))); };
            auto eta_label = [&](int iz){ return std::string(Form("%.1f < #eta^{pair} < %.1f",
                HO.num->GetZaxis()->GetBinLowEdge(iz), HO.num->GetZaxis()->GetBinUpEdge(iz))); };

            if (R.tag == "full") stat_full.assign(npt, std::vector<CellStat>(neta));

            // subplot grid: nrows >= ncols, nrows ~ sqrt(neta) (feedback_subplot_layout)
            const int ncol = (int)std::ceil(std::sqrt((double)neta));
            const int nrow = (int)std::ceil((double)neta / ncol);

            for (int iy = 1; iy <= npt; ++iy) {            // ONE PNG PER PAIR-pT BIN
                TCanvas c(Form("c_%s_%s_pt%d", S.key.c_str(), R.tag.c_str(), iy), "",
                          620 * ncol, 580 * nrow);
                c.Divide(ncol, nrow);
                for (int iz = 1; iz <= neta; ++iz) {       // ONE SUBPLOT PER PAIR-eta BIN
                    c.cd(iz);
                    auto pads = SplitPadForRatio(Form("%s_%s_%d_%d", S.key.c_str(),
                                                      R.tag.c_str(), iy, iz), false, 0.34);

                    TH1D* ro = CellRatio(HO.num, HO.den, HO.errA, HO.errB, HO.covP, HO.covQ, iy, iz,
                                         Form("o_%s_%s_%s_%d_%d", sample.c_str(), S.key.c_str(),
                                              R.tag.c_str(), iy, iz));
                    TH1D* rc = CellRatio(HC.num, HC.den, HC.errA, HC.errB, HC.covP, HC.covQ, iy, iz,
                                         Form("c_%s_%s_%s_%d_%d", sample.c_str(), S.key.c_str(),
                                              R.tag.c_str(), iy, iz));

                    double ymax = 0.;
                    for (int i = 1; i <= ro->GetNbinsX(); ++i) {
                        ymax = std::max(ymax, ro->GetBinContent(i) + ro->GetBinError(i));
                        ymax = std::max(ymax, rc->GetBinContent(i) + rc->GetBinError(i));
                    }
                    // CAP the auto-range (a single noisy near-empty cell would squeeze the
                    // structure this panel exists to show); every off-scale point is ARROWED and
                    // listed on the canvas -- silently dropping data is not allowed.
                    ymax = std::min(std::max(1.15, 1.15 * ymax), 3.0);

                    pads.first->cd();
                    DrawEffFrame(0.0, R.xhi, "", 0.0, ymax, S.ytitle, true);
                    DrawUnityLine(0.0, R.xhi);
                    ro->SetMarkerStyle(20); ro->SetMarkerColor(kOrigColor);
                    ro->SetLineColor(kOrigColor); ro->SetLineWidth(2);
                    rc->SetMarkerStyle(21); rc->SetMarkerColor(kCorrColor);
                    rc->SetLineColor(kCorrColor); rc->SetLineWidth(2);
                    ro->Draw("E1 same");
                    rc->Draw("E1 same");
                    auto offo = MarkOffScaleHist(ro, 0.0, ymax, kOrigColor);
                    auto offc = MarkOffScaleHist(rc, 0.0, ymax, kCorrColor);
                    offo.insert(offo.end(), offc.begin(), offc.end());

                    TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextSize(0.055);
                    tl.DrawLatex(0.17, 0.88, eta_label(iz).c_str());
                    if (!offo.empty()) {
                        TLatex nt; nt.SetNDC(); nt.SetTextFont(42); nt.SetTextSize(0.032);
                        nt.SetTextColor(kGray + 3);
                        double y = 0.24;
                        for (size_t i = 0; i < offo.size() && i < 6; ++i) {
                            nt.DrawLatex(0.17, y, ((i == 0 ? std::string("off scale (arrows): ")
                                                          : std::string("  ")) + offo[i]).c_str());
                            y -= 0.045;
                        }
                        if (offo.size() > 6)
                            nt.DrawLatex(0.17, y, Form("  ... and %d more", (int)offo.size() - 6));
                    }
                    if (iz == 1) {
                        auto* leg = new TLegend(0.40, 0.62, 0.95, 0.86);
                        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.045);
                        leg->AddEntry(ro, "original MC", "lp");
                        leg->AddEntry(rc, "corrected MC", "lp");
                        leg->Draw();
                    }

                    // ---- PULL pad: (corrected - original) / sigma_original --------------------
                    // A RATIO pad would need an error, and the two series come from the SAME
                    // events with the same trigger decisions -- their difference is systematic,
                    // not statistical, so a quadrature error would be meaningless. The pull
                    // against the ORIGINAL measurement's own conditional error is the quantity
                    // that actually answers Q2: "is the shift small compared with the precision
                    // of the correction being shifted?"
                    pads.second->cd();
                    auto* pull = (TH1D*)ro->Clone(Form("pull_%s_%s_%s_%d_%d", sample.c_str(),
                                                       S.key.c_str(), R.tag.c_str(), iy, iz));
                    pull->SetDirectory(nullptr);
                    double maxpull = 0., maxabs = 0.;
                    std::vector<double> pulls;
                    for (int i = 1; i <= pull->GetNbinsX(); ++i) {
                        const double o = ro->GetBinContent(i), cc = rc->GetBinContent(i);
                        const double e = ro->GetBinError(i);
                        if (e <= 0. || (o == 0. && cc == 0.)) { pull->SetBinContent(i, 0.);
                                                                pull->SetBinError(i, 0.); continue; }
                        const double p = (cc - o) / e;
                        pull->SetBinContent(i, p);
                        pull->SetBinError(i, 0.);
                        maxpull = std::max(maxpull, std::fabs(p));
                        maxabs  = std::max(maxabs, std::fabs(cc - o));
                        pulls.push_back(std::fabs(p));
                    }
                    double ylim = std::max(1.0, 1.25 * maxpull);
                    ylim = std::min(ylim, 5.0);
                    DrawRatioFrame(0.0, R.xhi, "#DeltaR",
                                   "#frac{corr - orig}{#sigma_{orig}}", -ylim, ylim, 0.0, 0.34);
                    pull->SetMarkerStyle(20);
                    pull->SetMarkerColor(kGreen + 2);
                    pull->SetLineColor(kGreen + 2);
                    pull->SetMarkerSize(0.8);
                    pull->Draw("P same");
                    // 2-specifier format: MarkOffScaleHist always passes (x, value, error), and
                    // the trailing error argument is simply unused here (a pull carries no error).
                    MarkOffScaleHist(pull, -ylim, ylim, kGreen + 2, "#DeltaR=%.2f: %.1f#sigma");

                    if (R.tag == "full") {
                        std::sort(pulls.begin(), pulls.end());
                        CellStat cs;
                        cs.maxpull = maxpull;
                        cs.medpull = pulls.empty() ? -1. : pulls[pulls.size() / 2];
                        cs.maxabs  = maxabs;
                        cs.nb      = (int)pulls.size();
                        cs.plat_o  = PlateauWeightedMean(ro, kPlateauLo, kPlateauHi).first;
                        cs.plat_c  = PlateauWeightedMean(rc, kPlateauLo, kPlateauHi).first;
                        stat_full[iy - 1][iz - 1] = cs;
                    }
                    c.cd(iz);
                }
                c.cd(0);
                TLatex st; st.SetNDC(); st.SetTextFont(42); st.SetTextSize(0.014);
                st.DrawLatex(0.02, 0.982, (headline + "  --  " + S.key + ", " + R.tag +
                             " #DeltaR, " + pt_label(iy) +
                             "  --  original vs SF-corrected MC, one panel per #eta^{pair} bin").c_str());
                SaveCanvas(c, S.dir + S.key + "_corr_vs_orig_" + R.tag + "_pairpt" +
                               std::to_string(iy) + ".png");
            }
        }

        // ---- per-cell pull table (full-dR range: it contains the plateau) ----
        {
            Step34Hists HO = LoadStep34(fo, S.base, "full", S.with_cov);
            const int npt  = (int)stat_full.size();
            const int neta = npt ? (int)stat_full[0].size() : 0;
            auto pt_hdr  = [&](int iy){ return Form("pTpair[%.0f,%.0f)",
                HO.num->GetYaxis()->GetBinLowEdge(iy), HO.num->GetYaxis()->GetBinUpEdge(iy)); };
            auto eta_hdr = [&](int iz){ return Form("[%.1f,%.1f)",
                HO.num->GetZaxis()->GetBinLowEdge(iz), HO.num->GetZaxis()->GetBinUpEdge(iz)); };

            std::ofstream os(S.dir + S.key + "_corr_vs_orig_pulls.txt");
            os << "# Q2: does correcting the MC single-muon efficiency to data change the "
               << S.key << " dR correction?\n";
            os << "# sample=" << sample << "  WP=" << wp_text << "  " << cfg.sample_text << "\n";
            os << "# Numerator weight: original w/eps_MC  vs  corrected w*SF/eps_corr.\n";
            os << "# If eps_corr were exactly eps_data, SF/eps_corr = 1/eps_MC and the two would be\n"
                  "# IDENTICAL bin by bin. Residuals can only come from (i) the corrected fit not\n"
                  "# being exactly the data fit, (ii) the eps floor/cap guards, (iii) the q.eta-gap\n"
                  "# 2D fallback.\n";
            os << "# Cells: max|pull| / median|pull| / max|corr-orig| , pull = (corr-orig)/sigma_orig\n";
            os << "# (sigma_orig = the round-7 CONDITIONAL error of the ORIGINAL measurement; the\n"
                  "#  two series share the same events, so a quadrature error would be meaningless.)\n";
            os << "# dR range: FULL (0-5.75), all " << (npt ? stat_full[0][0].nb : 0)
               << "-ish measured bins per cell.\n\n";
            os << std::left << std::setw(20) << "pair-eta \\ pair-pT";
            for (int iy = 1; iy <= npt; ++iy) os << std::setw(30) << pt_hdr(iy);
            os << "\n";
            double gmax = 0., gmedsum = 0.; int gn = 0;
            for (int iz = 1; iz <= neta; ++iz) {
                os << std::left << std::setw(20) << eta_hdr(iz);
                for (int iy = 1; iy <= npt; ++iy) {
                    const CellStat& s = stat_full[iy - 1][iz - 1];
                    if (s.nb <= 0) { os << std::setw(30) << "--"; continue; }
                    os << std::setw(30) << Form("%.3f/%.3f/%.4f", s.maxpull, s.medpull, s.maxabs);
                    gmax = std::max(gmax, s.maxpull);
                    gmedsum += s.medpull; ++gn;
                }
                os << "\n";
            }
            os << "\n# Large-dR plateau (weighted mean over dR in [" << kPlateauLo << ","
               << kPlateauHi << "]) per cell: original -> corrected\n";
            os << std::left << std::setw(20) << "pair-eta \\ pair-pT";
            for (int iy = 1; iy <= npt; ++iy) os << std::setw(30) << pt_hdr(iy);
            os << "\n";
            for (int iz = 1; iz <= neta; ++iz) {
                os << std::left << std::setw(20) << eta_hdr(iz);
                for (int iy = 1; iy <= npt; ++iy) {
                    const CellStat& s = stat_full[iy - 1][iz - 1];
                    os << std::setw(30) << ((s.plat_o > 0 && s.plat_c > 0)
                        ? Form("%.4f -> %.4f", s.plat_o, s.plat_c) : "--");
                }
                os << "\n";
            }
            os << "\n# VERDICT (" << S.key << "): max |pull| over all cells = "
               << Form("%.3f", gmax) << " sigma ; mean of the per-cell median |pull| = "
               << Form("%.3f", gn ? gmedsum / gn : -1.) << " sigma\n";
            std::cout << "  " << S.key << " Q2: max |pull| = " << gmax
                      << " sigma, mean per-cell median |pull| = " << (gn ? gmedsum / gn : -1.)
                      << " sigma\n";
            std::cout << "  wrote " << S.dir << S.key << "_corr_vs_orig_pulls.txt\n";
        }

        fo->Close();
        fc->Close();
    }

    fcorr1->Close(); fnom1->Close(); fcfit->Close(); fnfit->Close();
    fdata->Close();  fdfit->Close();
    std::cout << "\nplot_mc_trig_eff_corrected(" << sample << ", "
              << (use_tight_wp ? "tight" : "medium") << ") done.\n";
}
