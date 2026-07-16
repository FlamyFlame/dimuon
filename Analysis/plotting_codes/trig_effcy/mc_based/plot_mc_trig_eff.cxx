// plot_mc_trig_eff.cxx
// Step-1/2/3 MC trigger-efficiency plots (mc_trigger_efficiency.md §3.1-§3.3).
//
//   Step 1 (§3.1): MC direct P[mu4-matched | offline reco muon] vs data tag-and-probe
//                  P[2mu4 | mu4-tag, DeltaR>0.8] (the data eps^nc estimate), per charge.
//   Step 2 (§3.2): DeltaR-binned MC singles efficiency (factorization cross-check),
//                  3 DeltaR bins overlaid + Step-1 inclusive reference.
//   Step 3 (§3.3): eps_dR(dR) = inverse-weighted num / unweighted denom (TH1::Divide with
//                  error propagation -- weights > 1 make Bayes invalid, same convention as
//                  the data cross-term), plateau = weighted mean over dR in [1,4].
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

#include <cctype>
#include <cmath>
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
    } else {
        throw std::runtime_error("plot_mc_trig_eff: sample must be 'pp' or 'overlay', got " + sample);
    }
    return c;
}

// charge <-> data-sign mapping (F3): sign1 = mu+, sign2 = mu-
const std::vector<std::string> kCharges     = {"muplus", "muminus"};
const std::vector<std::string> kDataSigns   = {"sign1", "sign2"};
const std::vector<std::string> kChargeTex   = {"#mu^{+}", "#mu^{-}"};

// fine q.eta bins (CommonEffcyConfig.h q_eta_proj_ranges_fine_excl_gap)
const std::vector<std::string> kQEtaSuffix = {
    "minus2_40_TO_minus2_00", "minus2_00_TO_minus1_60", "minus1_60_TO_minus1_30",
    "minus0_90_TO_minus0_50", "minus0_50_TO_minus0_10", "0_10_TO_0_50",
    "0_50_TO_1_00", "1_30_TO_1_60", "1_60_TO_2_00", "2_00_TO_2_20"};
const std::vector<std::pair<double,double>> kQEtaRange = {
    {-2.4,-2.0},{-2.0,-1.6},{-1.6,-1.3},{-0.9,-0.5},{-0.5,-0.1},
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
    const std::string dir2 = cfg.out_base + "step2_dr_binned_singles/" + wp_dir;
    const std::string dir3 = cfg.out_base + "step3_dr_correction/"     + wp_dir;
    for (const auto& d : {dir1, dir2, dir3}) gSystem->mkdir(d.c_str(), kTRUE);

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
            auto* gda = GetObj<TGraphAsymmErrors>(fdata,
                "g_pt2nd_vs_q_eta2nd" + cfg.ctr + "_" + kDataSigns[ic] +
                "_2mu4_sepr_py_" + kQEtaSuffix[iq] + "_divided");

            gmc = (TGraphAsymmErrors*)gmc->Clone();
            gda = (TGraphAsymmErrors*)gda->Clone();
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
            gmc->Draw("PZ same");

            TLatex tl;
            tl.SetNDC();
            tl.SetTextSize(0.075);
            tl.SetTextFont(42);
            tl.DrawLatex(0.30, 0.12, Form("%.1f < q#upoint#eta < %.1f",
                                          kQEtaRange[iq].first, kQEtaRange[iq].second));

            // R3 ratio pad: MC / data per q.eta bin. This is where the (-2.4,-2.0) low-pT
            // discrepancy lives -- the ratio makes it quantitative instead of eyeballed.
            pads.second->cd();
            gPad->SetLogx();
            DrawRatioFrame(4.0, 60.0, "p_{T} [GeV]", "MC / data", 0.5, 2.6);
            auto* grat = DivideGraphClean(gmc, gda);
            StyleGraph(grat, kMCColor, 21, 0.7);
            grat->Draw("PZ same");
            MarkOffScale(grat, 0.5, 2.6, kMCColor);
            c.cd(static_cast<int>(iq) + 1);
        }
        // legend / label pad
        c.cd(11);
        DrawHeadline(headline + ", " + kChargeTex[ic], 0.02, 0.88, 0.048);
        auto* leg = new TLegend(0.02, 0.45, 0.98, 0.80);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.05);
        auto* gm = new TGraphAsymmErrors(); StyleGraph(gm, kMCColor, 21);
        auto* gd = new TGraphAsymmErrors(); StyleGraph(gd, kDataColor, 20);
        leg->AddEntry(gm, (mc_leg + " (+ fit)").c_str(), "lp");
        leg->AddEntry(gd, (cfg.data_text + " tag&probe").c_str(), "lp");
        leg->Draw();
        TLatex note;
        note.SetNDC();
        note.SetTextSize(0.045);
        note.SetTextFont(42);
        note.DrawLatex(0.02, 0.30, "Data: T&P P(2mu4 | mu4 tag, #DeltaR>0.8 pairs);");
        note.DrawLatex(0.02, 0.22, "MC: direct conditional, no T&P");
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
                rnum = GetObj<TH1D>(fmc, "h_mc_pt_num_"   + kCharges[ic]);
                rden = GetObj<TH1D>(fmc, "h_mc_pt_denom_" + kCharges[ic]);
            } else {
                TH2D* h2n = GetObj<TH2D>(fmc, "h_mc_pt_vs_q_eta_num_"   + kCharges[ic]);
                TH2D* h2d = GetObj<TH2D>(fmc, "h_mc_pt_vs_q_eta_denom_" + kCharges[ic]);
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
                TH1* n = GetObj<TH1D>(fmc, "h_mc_pair_" + v.tag + "_num_"   +
                                            kCharges[ic] + "_" + kDrSuffix[id]);
                TH1* d = GetObj<TH1D>(fmc, "h_mc_pair_" + v.tag + "_denom_" +
                                            kCharges[ic] + "_" + kDrSuffix[id]);
                n = Coarsen(n, Form("rb_n_%s_%d_%zu", v.tag.c_str(), ic, id));
                d = Coarsen(d, Form("rb_d_%s_%d_%zu", v.tag.c_str(), ic, id));
                auto* g = BayesEff(n, d);
                StyleGraph(g, kDrColor[id], kDrMarker[id], 0.8);
                leg->AddEntry(g, kDrTex[id].c_str(), "lp");
                gdr.push_back(g);
            }
            leg->AddEntry(gref, leg_incl.c_str(), "l");
            leg->Draw();                       // legend first ...
            gref->Draw("LX same");             // ... then the data on top of it
            for (auto* g : gdr) g->Draw("PZ same");
            DrawHeadline(headline + ", " + kChargeTex[ic], 0.14, 0.955, 0.05);

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
                TH2D* h2n = GetObj<TH2D>(fmc, "h_mc_pair_pt_vs_q_eta_num_" +
                                              kCharges[ic] + "_" + kDrSuffix[id]);
                TH2D* h2d = GetObj<TH2D>(fmc, "h_mc_pair_pt_vs_q_eta_denom_" +
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
        c.cd(11);
        DrawHeadline(headline + ", " + kChargeTex[ic], 0.02, 0.88, 0.048); // 0.048: long overlay headline + charge must fit (review iter 1)
        auto* leg = new TLegend(0.05, 0.35, 0.95, 0.80);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.06);
        for (size_t id = 0; id < kDrSuffix.size(); ++id) {
            auto* g = new TGraphAsymmErrors();
            StyleGraph(g, kDrColor[id], kDrMarker[id]);
            leg->AddEntry(g, kDrTex[id].c_str(), "lp");
        }
        leg->Draw();
        SaveCanvas(c, dir2 + "step2_eff_pt_in_q_eta_bins_" + kCharges[ic] +
                          "_dr_bins.png");
    }

    // ================================================================
    // Step 3 (§3.3): eps_dR(dR) = inverse-weighted num / denom
    // ================================================================
    std::cout << "\n===== Step 3 (" << sample << ", " << wp_text << ") =====\n";

    // ratios with TH1::Divide error propagation (inverse weights > 1 -> Bayes invalid)
    auto MakeRatio = [&](const std::string& tag) -> TH1D* {
        TH1D* num = GetObj<TH1D>(fmc3, "h_mc_dr_" + tag + "_num");
        TH1D* den = GetObj<TH1D>(fmc3, "h_mc_dr_" + tag + "_denom");
        auto* r = (TH1D*)num->Clone(("r_dr_" + tag).c_str());
        r->SetDirectory(nullptr);
        r->Divide(den);
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
            auto* r = (TH1D*)n->Clone(Form("s3_r_%zu", is));
            r->SetDirectory(nullptr);
            r->Divide(d);
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

    fmc->Close(); fmc3->Close(); ffit->Close(); fdata->Close();
    std::cout << "\nplot_mc_trig_eff(" << sample << ") done.\n";
}
