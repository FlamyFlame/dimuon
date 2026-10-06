// =================================================================================================
// plot_pp24_fullsim_pair_reco_eff_views.cxx
//
// Projection views of the pp24 fullsim PAIR reconstruction efficiency on the CANONICAL coarse
// cells, read straight from the RDF hist-filling output (not from the applied product), so that
// any truth pair category can be drawn. docs/tracking/pair_reco_eff_origin_test_and_views.md.
//
// eps_reco(cell) = N[truth pair in T, both reco-matched, pair passes the WP, RECO pair passes the
//                    data signal selection incl. the detector-gap cuts]  /  N[truth pair in T],
// T = truth OS pair in the truth signal region (mass window, pair pT; NO gap cut), binned in TRUTH
// (pair pT coarse, pair eta coarse, Z) with Z = dR (applied) or the truth pair mass (alternative).
//
// TWO PARTS
//  1. ORIGIN TEST (Z = dR only): all truth OS pairs (`_op`) vs from_same_b OS pairs (`_single_b`,
//     the applied sample), overlaid. `_single_b` is a SUBSET of `_op`, so the overlap of the two
//     curves is not a statistical test; the companion table compares `_single_b` with its
//     INDEPENDENT complement (`_op` minus `_single_b`) bin by bin and over all 3D cells.
//  2. VIEWS of one category: eps vs eta^pair in pT panels; vs pT in eta^pair panels; vs Z in pT
//     panels and vs pT in Z panels, each for all |eta^pair| and per canonical |eta^pair| group.
//
// Every curve is a ratio of PROJECTED numerator and denominator (never a projection of a ratio).
// Axes are verified against ParamsSet / CommonEffcyConfig; the |eta^pair| groups come from the
// canonical fold MakeDrEtaGroups(axis, true) -- nothing is retyped.
//
// Usage (from Analysis/plotting_codes/reco_effcy/):
//   root -l -b -q 'plot_pp24_fullsim_pair_reco_eff_views.cxx+(true, "origin_test")'
//   root -l -b -q 'plot_pp24_fullsim_pair_reco_eff_views.cxx+(true, "views", "_op", "dr")'
//   arguments: tight_WP, part ("origin_test" | "views"), category ("_op" | "_single_b"),
//              z ("dr" | "minv"), use_full_sample
// =================================================================================================

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TNamed.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../../Utilities/PairRecoEffDefinition.h"
#include "../trig_effcy/mc_based/dr_correction_sample_cfg.h"
#include "../trig_effcy/mc_based/dr_correction_cell_groups.h"   // MakeDrEtaGroups: the |eta| fold

namespace {

using Sel = std::vector<std::pair<int, int>>;   // union of inclusive bin ranges on one axis

const double kNoMeasure = -999.0;
const std::vector<Color_t> kCols = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};
const std::vector<Style_t> kMrks = {20, 24, 21, 22, 23};

// The third axis of the 3D: dR (applied) or the truth pair mass (alternative).
struct ZVar {
    std::string key, hist_tag, sym, unit, ascii;
};
ZVar GetZVar(const std::string& z)
{
    if (z == "dr")   return {"dr", "truth_dr_effcy", "#DeltaR", "", "dR"};
    if (z == "minv") return {"minv", "truth_minv_effcy", "m_{#mu#mu}", " GeV", "m_mumu"};
    throw std::runtime_error("plot_pp24_fullsim_pair_reco_eff_views: unknown z '" + z + "'");
}

std::string HistName(const ZVar& zv, const std::string& cat, const std::string& leg)
{
    return "h_" + zv.hist_tag + "_vs_truth_pair_eta_coarse_vs_truth_pair_pt_coarse" + cat + leg;
}

TH3D* Get3D(TFile* f, const std::string& n)
{
    auto* h = dynamic_cast<TH3D*>(f->Get(n.c_str()));
    if (!h) throw std::runtime_error("plot_pp24_fullsim_pair_reco_eff_views: missing " + n
                                     + " in " + f->GetName());
    return h;
}

void CheckAxis(const TAxis* ax, const std::vector<double>& edges, const std::string& what)
{
    if (ax->GetNbins() != (int)edges.size() - 1)
        throw std::runtime_error("plot_pp24_fullsim_pair_reco_eff_views: " + what + " axis has "
                                 + std::to_string(ax->GetNbins()) + " bins, canonical "
                                 + std::to_string(edges.size() - 1));
    for (int i = 1; i <= ax->GetNbins() + 1; ++i)
        if (std::fabs(ax->GetBinLowEdge(i) - edges[i - 1]) > 1e-5)
            throw std::runtime_error("plot_pp24_fullsim_pair_reco_eff_views: " + what + " edge "
                                     + std::to_string(i - 1) + " differs from the canonical one");
}

void CheckCanonical(const TH3D* h, const ZVar& zv)
{
    static const ParamsSet pms{};
    static const CommonEffcyConfig cfg{};
    std::vector<double> eta{cfg.pair_eta_proj_ranges_coarse_incl_gap.front().first};
    for (const auto& r : cfg.pair_eta_proj_ranges_coarse_incl_gap) eta.push_back(r.second);
    CheckAxis(h->GetXaxis(), pms.pair_pt_coarse_bins, "pair-pT");
    CheckAxis(h->GetYaxis(), eta, "pair-eta");
    if (zv.key == "dr") {
        std::vector<double> dr(cfg.dr_bins_edges_for_reco_effcy.begin(),
                               cfg.dr_bins_edges_for_reco_effcy.end());
        CheckAxis(h->GetZaxis(), dr, "dR");
    }
}

// Complement of a subset: content difference, Sum(w^2) difference (the subset's entries are a
// subset of the superset's, so their variances subtract -- TH1::Add(h, -1) would ADD them).
TH3D* Complement(const TH3D* sup, const TH3D* sub, const std::string& name)
{
    auto* c = static_cast<TH3D*>(sup->Clone(name.c_str()));
    c->SetDirectory(nullptr);
    for (int b = 0; b < sup->GetNcells(); ++b) {
        const double v  = sup->GetBinContent(b) - sub->GetBinContent(b);
        const double e2 = std::pow(sup->GetBinError(b), 2) - std::pow(sub->GetBinError(b), 2);
        c->SetBinContent(b, v);
        c->SetBinError(b, std::sqrt(std::max(0.0, e2)));
    }
    return c;
}

Sel All(int n) { return {{1, n}}; }
Sel One(int i) { return {{i, i}}; }

// 1D projection of `h` onto axis `a` ('x','y','z'), summing the OTHER two axes over the unions
// s1, s2 (in axis order, skipping `a`). Sum(w^2) carried exactly.
TH1D* Project(const TH3D* h, char a, const Sel& s1, const Sel& s2, const std::string& name)
{
    const TAxis* ax = (a == 'x') ? h->GetXaxis() : (a == 'y') ? h->GetYaxis() : h->GetZaxis();
    std::vector<double> edges;
    for (int i = 1; i <= ax->GetNbins() + 1; ++i) edges.push_back(ax->GetBinLowEdge(i));
    auto* p = new TH1D(name.c_str(), "", ax->GetNbins(), edges.data());
    p->SetDirectory(nullptr);
    p->Sumw2();
    auto in = [](int i, const Sel& s) {
        for (const auto& r : s) if (i >= r.first && i <= r.second) return true;
        return false;
    };
    for (int i = 1; i <= h->GetNbinsX(); ++i)
        for (int j = 1; j <= h->GetNbinsY(); ++j)
            for (int k = 1; k <= h->GetNbinsZ(); ++k) {
                int keep, o1, o2;
                if (a == 'x')      { keep = i; o1 = j; o2 = k; }
                else if (a == 'y') { keep = j; o1 = i; o2 = k; }
                else               { keep = k; o1 = i; o2 = j; }
                if (!in(o1, s1) || !in(o2, s2)) continue;
                const int b = h->GetBin(i, j, k);
                p->SetBinContent(keep, p->GetBinContent(keep) + h->GetBinContent(b));
                p->SetBinError(keep, std::hypot(p->GetBinError(keep), h->GetBinError(b)));
            }
    return p;
}

// eps = num/den with the binomial error for weighted histograms (TH1::Divide "B"). Empty
// denominator -> no measure (pushed off the frame). eps = 0 or 1 with zero error -> drawn with the
// 1-sigma Wilson width 1/(Neff+1), as in plot_pp24_fullsim_pair_reco_eff.cxx.
TH1D* Eff(const TH3D* num, const TH3D* den, char a, const Sel& s1, const Sel& s2,
          const std::string& name)
{
    TH1D* n = Project(num, a, s1, s2, name + "_n");
    TH1D* d = Project(den, a, s1, s2, name + "_d");
    auto* e = static_cast<TH1D*>(n->Clone(name.c_str()));
    e->SetDirectory(nullptr);
    e->Divide(n, d, 1.0, 1.0, "B");
    for (int b = 1; b <= e->GetNbinsX(); ++b) {
        const double dc = d->GetBinContent(b), de = d->GetBinError(b);
        if (dc <= 0.) { e->SetBinContent(b, kNoMeasure); e->SetBinError(b, 0.); continue; }
        const double v = e->GetBinContent(b);
        if (e->GetBinError(b) <= 0. && !(v > 0. && v < 1.) && de > 0.)
            e->SetBinError(b, 1.0 / (std::pow(dc / de, 2) + 1.0));
    }
    delete n;
    delete d;
    return e;
}

bool Measured(const TH1* h, int b) { return h->GetBinContent(b) > kNoMeasure / 2.0; }

std::string RangeLabel(const std::string& sym, const std::string& unit, double lo, double hi,
                       int prec)
{
    return Form("%.*f < %s < %.*f%s", prec, lo, sym.c_str(), prec, hi, unit.c_str());
}
std::string PtLabel(const TAxis* a, int i)
{
    return Form("%.4g < p_{T}^{pair} < %.4g GeV", a->GetBinLowEdge(i), a->GetBinUpEdge(i));
}
std::string EtaLabel(const TAxis* a, int j)
{
    return RangeLabel("#eta^{pair}", "", a->GetBinLowEdge(j), a->GetBinUpEdge(j), 1);
}
std::string ZLabel(const ZVar& zv, const TAxis* a, int k)
{
    return RangeLabel(zv.sym, zv.unit, a->GetBinLowEdge(k), a->GetBinUpEdge(k),
                      zv.key == "dr" ? 1 : 2);
}
std::string AbsEtaLabel(const DrAxisGroups& G, int g)
{
    const double lo = G.edges[g], hi = G.edges[g + 1];
    if (lo <= 0.) return Form("|#eta^{pair}| < %.1f", hi);
    return Form("%.1f < |#eta^{pair}| < %.1f", lo, hi);
}
std::string AbsEtaTag(const DrAxisGroups& G, int g)
{
    auto f = [](double x) { std::string s = Form("%.1f", x); std::replace(s.begin(), s.end(), '.', 'p'); return s; };
    return "abs_pair_eta_" + f(G.edges[g]) + "_" + f(G.edges[g + 1]);
}
std::string Ascii(std::string s)
{
    for (const char* from : {"#eta^{pair}", "p_{T}^{pair}", "#DeltaR", "m_{#mu#mu}"}) {
        const std::string to = (std::string(from) == "#eta^{pair}") ? "eta_pair"
                             : (std::string(from) == "p_{T}^{pair}") ? "pT_pair"
                             : (std::string(from) == "#DeltaR") ? "dR" : "m_mumu";
        for (size_t p; (p = s.find(from)) != std::string::npos;) s.replace(p, std::string(from).size(), to);
    }
    return s;
}

// ------------------------------------------------------------------------------------------------
// One multi-panel figure: panel p holds curves[p] (one per series), x title, log x.
// One legend for the whole figure in the header strip (a per-panel legend covers points).
struct Series { std::string label; };
struct Figure {
    std::string path, xtitle, header, definition, wp_text, extra_header;
    bool logx = false;
    std::vector<std::string> panel_labels;
    std::vector<std::vector<TH1D*>> curves;    // [panel][series]
    std::vector<Series> series;
    std::vector<TH1D*> ratio;                  // optional, one per panel (bottom pad)
    std::string ratio_title;
};

// r = eps_sb / eps_op with eps_op = f*eps_sb + (1-f)*eps_c (f = single-b share of the projected
// denominator). The two overlaid curves share 85 % of their pairs, so their naive ratio error is
// wrong; the error is propagated through the INDEPENDENT pieces eps_sb and eps_c instead.
TH1D* OriginRatio(const TH1D* sb, const TH1D* co, const TH1D* d_sb, const TH1D* d_op,
                  const std::string& name)
{
    auto* r = static_cast<TH1D*>(sb->Clone(name.c_str()));
    r->SetDirectory(nullptr);
    for (int b = 1; b <= r->GetNbinsX(); ++b) {
        r->SetBinContent(b, kNoMeasure); r->SetBinError(b, 0.);
        if (!Measured(sb, b) || d_op->GetBinContent(b) <= 0.) continue;
        const double f  = d_sb->GetBinContent(b) / d_op->GetBinContent(b);
        const double es = sb->GetBinContent(b), ss = sb->GetBinError(b);
        const double ec = Measured(co, b) ? co->GetBinContent(b) : 0., sc = Measured(co, b) ? co->GetBinError(b) : 0.;
        const double eo = f * es + (1. - f) * ec;
        if (eo <= 0.) continue;
        const double ds = 1. / eo - es * f / (eo * eo), dc = -es * (1. - f) / (eo * eo);
        r->SetBinContent(b, es / eo);
        r->SetBinError(b, std::sqrt(ds * ds * ss * ss + dc * dc * sc * sc));
    }
    return r;
}

void DrawFigure(const Figure& F)
{
    const int np = (int)F.curves.size();
    // nrows >= ncols, nrows ~ sqrt(N); N <= 3 -> one row
    const int ncol = (np <= 3) ? np : (int)std::ceil(std::sqrt((double)np));
    const int nrow = (int)std::ceil((double)np / ncol);
    const int hstrip = F.extra_header.empty() ? 115 : 150;
    TCanvas c("c_fig", "", 430 * ncol, 340 * nrow + hstrip);
    const double fh = (double)hstrip / (340.0 * nrow + hstrip);
    TPad* head = new TPad("head", "", 0.0, 1.0 - fh, 1.0, 1.0);
    TPad* body = new TPad("body", "", 0.0, 0.0, 1.0, 1.0 - fh);
    head->SetFillStyle(0);
    body->SetFillStyle(0);
    head->Draw();
    body->Draw();
    head->cd();
    {
        const double sz = (ncol <= 2 ? 0.20 : 0.24) * 110.0 / hstrip;
        TLatex t;
        t.SetNDC();
        t.SetTextFont(42);
        t.SetTextSize(sz);
        double y = 0.80;
        t.DrawLatex(0.02, y, F.header.c_str());
        t.SetTextAlign(31);
        t.DrawLatex(0.98, y, F.wp_text.c_str());
        t.SetTextAlign(11);
        y -= 1.35 * sz;
        t.DrawLatex(0.02, y, F.definition.c_str());
        if (!F.extra_header.empty()) { y -= 1.35 * sz; t.DrawLatex(0.02, y, F.extra_header.c_str()); }
        TLegend* leg = new TLegend(0.02, 0.02, 0.98, y - 0.6 * sz);
        leg->SetNColumns((int)F.series.size());
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(sz);
        for (size_t s = 0; s < F.series.size(); ++s)
            for (int p = 0; p < np; ++p)
                if (F.curves[p][s]) { leg->AddEntry(F.curves[p][s], F.series[s].label.c_str(), "lp"); break; }
        leg->Draw();
    }
    body->cd();
    body->Divide(ncol, nrow);
    for (int p = 0; p < np; ++p) {
        body->cd(p + 1);
        TVirtualPad* cell = gPad;
        const bool with_ratio = !F.ratio.empty() && F.ratio[p];
        TPad* top = nullptr;
        if (with_ratio) {
            top = new TPad(Form("top%d", p), "", 0.0, 0.32, 1.0, 1.0);
            top->SetBottomMargin(0.02);
            top->Draw();
            top->cd();
        }
        gPad->SetTicks(1, 1);
        gPad->SetLogx(F.logx);
        gPad->SetLeftMargin(0.14);
        if (!with_ratio) gPad->SetBottomMargin(0.14);
        gPad->SetTopMargin(with_ratio ? 0.15 : 0.12);
        bool first = true;
        for (size_t s = 0; s < F.curves[p].size(); ++s) {
            TH1D* h = F.curves[p][s];
            if (!h) continue;
            h->SetLineColor(kCols[s % kCols.size()]);
            h->SetMarkerColor(kCols[s % kCols.size()]);
            h->SetMarkerStyle(kMrks[s % kMrks.size()]);
            h->SetMarkerSize(1.1);
            h->SetLineWidth(2);
            h->SetStats(0);
            h->SetTitle("");
            h->GetXaxis()->SetTitle(F.xtitle.c_str());
            h->GetXaxis()->SetNoExponent();
            h->GetXaxis()->SetMoreLogLabels(F.logx);
            h->GetYaxis()->SetTitle("#varepsilon_{reco}^{pair}");
            h->GetYaxis()->SetRangeUser(0.0, 1.1);
            h->GetXaxis()->SetTitleSize(0.052);
            h->GetXaxis()->SetLabelSize(0.044);
            h->GetYaxis()->SetTitleSize(0.052);
            h->GetYaxis()->SetLabelSize(0.044);
            h->GetYaxis()->SetTitleOffset(1.2);
            h->Draw(first ? "PE" : "PE SAME");
            first = false;
        }
        TLatex t;
        t.SetNDC();
        t.SetTextFont(42);
        t.SetTextSize(with_ratio ? 0.070 : 0.052);
        t.DrawLatex(0.17, with_ratio ? 0.89 : 0.915, F.panel_labels[p].c_str());   // above the frame, in the top margin
        if (with_ratio) {
            for (size_t s = 0; s < F.curves[p].size(); ++s)
                if (F.curves[p][s]) { F.curves[p][s]->GetXaxis()->SetLabelSize(0); F.curves[p][s]->GetXaxis()->SetTitleSize(0);
                                      F.curves[p][s]->GetYaxis()->SetTitleSize(0.070); F.curves[p][s]->GetYaxis()->SetLabelSize(0.060);
                                      F.curves[p][s]->GetYaxis()->SetTitleOffset(0.95); }
            cell->cd();
            TPad* bot = new TPad(Form("bot%d", p), "", 0.0, 0.0, 1.0, 0.32);
            bot->SetTopMargin(0.03);
            bot->SetBottomMargin(0.34);
            bot->SetLeftMargin(0.14);
            bot->SetTicks(1, 1);
            bot->SetLogx(F.logx);
            bot->Draw();
            bot->cd();
            TH1D* r = F.ratio[p];
            r->SetStats(0);
            r->SetTitle("");
            r->SetLineColor(kRed + 1); r->SetMarkerColor(kRed + 1); r->SetMarkerStyle(24); r->SetLineWidth(2);
            r->GetYaxis()->SetRangeUser(0.9, 1.1);
            r->GetYaxis()->SetNdivisions(504);
            r->GetYaxis()->SetTitle(F.ratio_title.c_str());
            r->GetYaxis()->SetTitleSize(0.13); r->GetYaxis()->SetLabelSize(0.12); r->GetYaxis()->SetTitleOffset(0.50);
            r->GetXaxis()->SetTitle(F.xtitle.c_str());
            r->GetXaxis()->SetNoExponent();
            r->GetXaxis()->SetMoreLogLabels(F.logx);
            r->GetXaxis()->SetTitleSize(0.14); r->GetXaxis()->SetLabelSize(0.12); r->GetXaxis()->SetTitleOffset(1.0);
            r->Draw("PE");
            TLine* one = new TLine(r->GetXaxis()->GetXmin(), 1.0, r->GetXaxis()->GetXmax(), 1.0);
            one->SetLineStyle(2);
            one->Draw();
        }
    }
    gSystem->mkdir(gSystem->DirName(F.path.c_str()), kTRUE);
    c.SaveAs(F.path.c_str());
}

// Values table: one block per figure, every drawn number.
void WriteTable(std::ofstream& tab, const Figure& F, const TAxis* xax)
{
    tab << "\n### " << gSystem->BaseName(F.path.c_str()) << "  (x = " << Ascii(F.xtitle) << ")\n";
    for (size_t p = 0; p < F.curves.size(); ++p) {
        tab << "  [" << Ascii(F.panel_labels[p]) << "]\n";
        for (size_t s = 0; s < F.series.size(); ++s) {
            const TH1D* h = F.curves[p][s];
            if (!h) continue;
            tab << "    " << Ascii(F.series[s].label) << ":";
            for (int b = 1; b <= h->GetNbinsX(); ++b)
                tab << "  [" << xax->GetBinLowEdge(b) << "," << xax->GetBinUpEdge(b) << "] "
                    << (Measured(h, b) ? std::string(Form("%.4f+-%.4f", h->GetBinContent(b), h->GetBinError(b)))
                                       : std::string("-"));
            tab << "\n";
        }
        if (!F.ratio.empty() && F.ratio[p]) {
            const TH1D* r = F.ratio[p];
            tab << "    ratio " << F.ratio_title << ":";
            for (int b = 1; b <= r->GetNbinsX(); ++b)
                tab << "  [" << xax->GetBinLowEdge(b) << "," << xax->GetBinUpEdge(b) << "] "
                    << (Measured(r, b) ? std::string(Form("%.4f+-%.4f", r->GetBinContent(b), r->GetBinError(b)))
                                       : std::string("-"));
            tab << "\n";
        }
    }
}

// Independent-sample comparison of two efficiency histograms: chi2 over bins measured in both.
struct Chi2 { double chi2 = 0; int ndf = 0; double max_pull = 0; };
void Accumulate(Chi2& c, const TH1D* a, const TH1D* b)
{
    for (int i = 1; i <= a->GetNbinsX(); ++i) {
        if (!Measured(a, i) || !Measured(b, i)) continue;
        const double s2 = std::pow(a->GetBinError(i), 2) + std::pow(b->GetBinError(i), 2);
        if (s2 <= 0.) continue;
        const double pull = (a->GetBinContent(i) - b->GetBinContent(i)) / std::sqrt(s2);
        c.chi2 += pull * pull;
        ++c.ndf;
        c.max_pull = std::max(c.max_pull, std::fabs(pull));
    }
}
// Same, restricted to bins whose denominator effective entry count N_eff = (sum w)^2 / sum w^2 is at
// least `min_neff` in BOTH samples: at small N_eff the Gaussian weighted-binomial error is not
// reliable (it shrinks towards 0 as eps -> 0 or 1) and a handful of such cells dominates the chi2.
double Neff(const TH1D* d, int b)
{
    return d->GetBinError(b) > 0. ? std::pow(d->GetBinContent(b) / d->GetBinError(b), 2) : 0.;
}
void AccumulateNeff(Chi2& c, const TH1D* a, const TH1D* b, const TH1D* da, const TH1D* db, double min_neff)
{
    for (int i = 1; i <= a->GetNbinsX(); ++i) {
        if (!Measured(a, i) || !Measured(b, i)) continue;
        if (Neff(da, i) < min_neff || Neff(db, i) < min_neff) continue;
        const double s2 = std::pow(a->GetBinError(i), 2) + std::pow(b->GetBinError(i), 2);
        if (s2 <= 0.) continue;
        const double pull = (a->GetBinContent(i) - b->GetBinContent(i)) / std::sqrt(s2);
        c.chi2 += pull * pull;
        ++c.ndf;
        c.max_pull = std::max(c.max_pull, std::fabs(pull));
    }
}
std::string Chi2Str(const Chi2& c)
{
    return Form("chi2/ndf = %.1f/%d, p = %.3f, max |pull| = %.2f", c.chi2, c.ndf,
                c.ndf ? TMath::Prob(c.chi2, c.ndf) : 1.0, c.max_pull);
}

}  // namespace

// =================================================================================================
void plot_pp24_fullsim_pair_reco_eff_views(bool tight_WP = true, std::string part = "origin_test",
                                           std::string category = "_op", std::string z = "dr",
                                           bool use_full_sample = true)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const ZVar zv = GetZVar(z);
    if (part == "origin_test" && zv.key != "dr")
        throw std::runtime_error("plot_pp24_fullsim_pair_reco_eff_views: the origin test is defined on the applied dR axis");
    if (category != "_op" && category != "_single_b")
        throw std::runtime_error("plot_pp24_fullsim_pair_reco_eff_views: category must be _op or _single_b (no SS reco-eff, user 2026-10-06)");

    const DrCorrSample id = GetDrCorrSample(use_full_sample ? "pp_full" : "pp", tight_WP);
    const std::string in_path = id.sample_dir + "histograms_pythia_fullsim_pp24_no_data_resonance_cuts"
                              + (use_full_sample ? "_full" : "") + ".root";
    const std::string wp      = tight_WP ? "tight" : "medium";
    const std::string wp_text = tight_WP ? "Tight muons" : "Medium muons";
    const std::string header  = std::string("Pythia8 full simulation, pp 2024 conditions")
                              + (use_full_sample ? "" : " (test sample)");
    const std::string out_root = id.sample_dir + "plots/pp24_reco_effcy_plots/" + wp + "/signal_region/";

    std::unique_ptr<TFile> fin(TFile::Open(in_path.c_str(), "READ"));
    if (!fin || fin->IsZombie()) throw std::runtime_error("cannot open " + in_path);
    {
        auto* def = dynamic_cast<TNamed*>(fin->Get(PairRecoEffDefinition::Key()));
        if (!def || std::string(def->GetTitle()) != PairRecoEffDefinition::Value())
            throw std::runtime_error(in_path + ": missing/different pair_reco_eff_definition marker -- "
                                     "rerun the fullsim RDF hist filling (stage 5)");
    }
    const std::string num_leg = "_pass_" + wp + "_and_signal_truth_and_reco";
    const std::string den_leg = "_pass_signal_truth";

    auto load = [&](const std::string& cat) {
        TH3D* n = Get3D(fin.get(), HistName(zv, cat, num_leg));
        TH3D* d = Get3D(fin.get(), HistName(zv, cat, den_leg));
        CheckCanonical(n, zv);
        CheckCanonical(d, zv);
        return std::make_pair(n, d);
    };

    // The pair category is named by the legend; the definition line stays short enough for a
    // two-column canvas.
    const auto def_text = [](const std::string&) {
        return std::string("#varepsilon_{reco}^{pair} = N_{reco} / N_{gen}; truth signal region; "
                           "N_{reco} incl. detector-gap cuts");
    };

    // ============================================================================================
    if (part == "origin_test") {
        auto [n_op, d_op] = load("_op");
        auto [n_sb, d_sb] = load("_single_b");
        TH3D* n_c = Complement(n_op, n_sb, "n_compl");
        TH3D* d_c = Complement(d_op, d_sb, "d_compl");
        const TAxis *apt = d_op->GetXaxis(), *aeta = d_op->GetYaxis(), *adr = d_op->GetZaxis();
        const int npt = apt->GetNbins(), neta = aeta->GetNbins(), ndr = adr->GetNbins();
        const std::string outdir = out_root + "origin_test/";
        gSystem->mkdir(outdir.c_str(), kTRUE);
        std::ofstream tab(outdir + "origin_test_values_" + wp + "WP.txt");
        tab << "Origin test of eps_reco (pp24 fullsim FULL, " << wp_text << "), input " << in_path << "\n"
            << "Series: all truth OS pairs in the signal region (_op) vs from_same_b OS pairs (_single_b).\n"
            << "_single_b is a SUBSET of _op: the statistical test compares _single_b with the INDEPENDENT\n"
            << "complement (_op - _single_b), pulls = (eps_sb - eps_compl)/sqrt(s_sb^2 + s_compl^2).\n";

        Figure FA, FB;
        for (Figure* F : {&FA, &FB}) {
            F->header = header;
            F->definition = def_text("");
            F->wp_text = wp_text;
            F->series = {{"all opposite-sign pairs"}, {"single-b opposite-sign pairs"}};
            F->ratio_title = "single-b / all";
        }
        FA.path = outdir + "pair_pt_in_pair_eta_bins_" + wp + "WP.png";
        FA.xtitle = "p_{T}^{pair} [GeV]";
        FA.logx = true;
        FA.extra_header = Form("#DeltaR < %.1f (integrated)", adr->GetBinUpEdge(ndr));
        FB.path = outdir + "dr_in_pair_pt_bins_" + wp + "WP.png";
        FB.xtitle = "#DeltaR";
        FB.extra_header = Form("%.1f < #eta^{pair} < %.1f (integrated)", aeta->GetBinLowEdge(1), aeta->GetBinUpEdge(neta));

        Chi2 cA, cB, c3;
        for (int j = 1; j <= neta; ++j) {
            TH1D* op = Eff(n_op, d_op, 'x', One(j), All(ndr), Form("A_op_%d", j));
            TH1D* sb = Eff(n_sb, d_sb, 'x', One(j), All(ndr), Form("A_sb_%d", j));
            TH1D* co = Eff(n_c, d_c, 'x', One(j), All(ndr), Form("A_co_%d", j));
            Chi2 cp; Accumulate(cp, sb, co); Accumulate(cA, sb, co);
            FA.panel_labels.push_back(EtaLabel(aeta, j));
            FA.curves.push_back({op, sb});
            {
                TH1D* dsb = Project(d_sb, 'x', One(j), All(ndr), Form("A_dsb_%d", j));
                TH1D* dop = Project(d_op, 'x', One(j), All(ndr), Form("A_dop_%d", j));
                FA.ratio.push_back(OriginRatio(sb, co, dsb, dop, Form("A_r_%d", j)));
                delete dsb; delete dop;
            }
            tab << "  [A] " << Ascii(EtaLabel(aeta, j)) << " : sb vs complement " << Chi2Str(cp) << "\n";
        }
        for (int i = 1; i <= npt; ++i) {
            TH1D* op = Eff(n_op, d_op, 'z', One(i), All(neta), Form("B_op_%d", i));
            TH1D* sb = Eff(n_sb, d_sb, 'z', One(i), All(neta), Form("B_sb_%d", i));
            TH1D* co = Eff(n_c, d_c, 'z', One(i), All(neta), Form("B_co_%d", i));
            Chi2 cp; Accumulate(cp, sb, co); Accumulate(cB, sb, co);
            FB.panel_labels.push_back(PtLabel(apt, i));
            FB.curves.push_back({op, sb});
            {
                TH1D* dsb = Project(d_sb, 'z', One(i), All(neta), Form("B_dsb_%d", i));
                TH1D* dop = Project(d_op, 'z', One(i), All(neta), Form("B_dop_%d", i));
                FB.ratio.push_back(OriginRatio(sb, co, dsb, dop, Form("B_r_%d", i)));
                delete dsb; delete dop;
            }
            tab << "  [B] " << Ascii(PtLabel(apt, i)) << " : sb vs complement " << Chi2Str(cp) << "\n";
        }
        // every 3D cell, the object actually applied; also with N_eff thresholds
        Chi2 c3_100, c3_1000;
        for (int j = 1; j <= neta; ++j)
            for (int k = 1; k <= ndr; ++k) {
                TH1D* sb = Eff(n_sb, d_sb, 'x', One(j), One(k), Form("C_sb_%d_%d", j, k));
                TH1D* co = Eff(n_c, d_c, 'x', One(j), One(k), Form("C_co_%d_%d", j, k));
                TH1D* dsb = Project(d_sb, 'x', One(j), One(k), Form("C_dsb_%d_%d", j, k));
                TH1D* dco = Project(d_c, 'x', One(j), One(k), Form("C_dco_%d_%d", j, k));
                Accumulate(c3, sb, co);
                AccumulateNeff(c3_100, sb, co, dsb, dco, 100.);
                AccumulateNeff(c3_1000, sb, co, dsb, dco, 1000.);
                delete sb; delete co; delete dsb; delete dco;
            }
        auto incl = [&](TH3D* n, TH3D* d) {
            double en, ed;
            const double N = n->IntegralAndError(1, npt, 1, neta, 1, ndr, en);
            const double D = d->IntegralAndError(1, npt, 1, neta, 1, ndr, ed);
            return std::make_pair(N / D, D * D / (ed * ed));
        };
        const auto [e_op, neff_op] = incl(n_op, d_op);
        const auto [e_sb, neff_sb] = incl(n_sb, d_sb);
        const auto [e_co, neff_co] = incl(n_c, d_c);
        tab << "\nSUMMARY (single-b vs independent complement)\n"
            << "  figure A (pT in eta^pair panels, dR-integrated): " << Chi2Str(cA) << "\n"
            << "  figure B (dR in pT panels, eta-integrated)     : " << Chi2Str(cB) << "\n"
            << "  all applied 3D cells                           : " << Chi2Str(c3) << "\n"
            << "  3D cells with N_eff >= 100 in both samples     : " << Chi2Str(c3_100) << "\n"
            << "  3D cells with N_eff >= 1000 in both samples    : " << Chi2Str(c3_1000) << "\n"
            << Form("  inclusive eps: _op %.4f (Neff %.0f), _single_b %.4f (Neff %.0f), complement %.4f (Neff %.0f)\n",
                    e_op, neff_op, e_sb, neff_sb, e_co, neff_co)
            << Form("  single-b weight fraction of the _op denominator: %.4f\n",
                    d_sb->Integral(1, npt, 1, neta, 1, ndr) / d_op->Integral(1, npt, 1, neta, 1, ndr));
        DrawFigure(FA);
        DrawFigure(FB);
        WriteTable(tab, FA, apt);
        WriteTable(tab, FB, adr);
        tab.close();
        std::cout << "[origin_test] " << wp << ": A " << Chi2Str(cA) << "\n                B " << Chi2Str(cB)
                  << "\n                3D " << Chi2Str(c3) << "\n                3D Neff>=100 " << Chi2Str(c3_100)
                  << "\n                3D Neff>=1000 " << Chi2Str(c3_1000) << "\n  -> " << outdir << std::endl;
        return;
    }

    // ============================================================================================
    if (part != "views") throw std::runtime_error("part must be origin_test or views");
    auto [num, den] = load(category);
    const TAxis *apt = den->GetXaxis(), *aeta = den->GetYaxis(), *azz = den->GetZaxis();
    const int npt = apt->GetNbins(), neta = aeta->GetNbins(), nz = azz->GetNbins();
    const DrAxisGroups G = MakeDrEtaGroups(aeta, true);   // canonical |eta^pair| fold
    const std::string sample_dir = (category == "_op") ? "all_os/" : "single_b/";
    const std::string outdir = out_root + sample_dir;
    gSystem->mkdir(outdir.c_str(), kTRUE);
    std::ofstream tab(outdir + "views_" + zv.key + "_values_" + wp + "WP.txt");
    tab << "eps_reco views (pp24 fullsim FULL, " << wp_text << ", category " << category
        << ", 3rd axis " << zv.ascii << "), input " << in_path << "\n";

    std::vector<Figure> figs;
    auto base = [&](const std::string& path, const std::string& xt, bool logx, const std::string& extra) {
        Figure F;
        F.path = path; F.xtitle = xt; F.logx = logx; F.extra_header = extra;
        F.header = header; F.definition = def_text(category); F.wp_text = wp_text;
        F.series = {{category == "_op" ? "all opposite-sign pairs" : "single-b opposite-sign pairs"}};
        return F;
    };
    const std::string z_int = Form("%s < %.*f%s (integrated)", zv.sym.c_str(), zv.key == "dr" ? 1 : 2,
                                   azz->GetBinUpEdge(nz), zv.unit.c_str());
    const std::string z_int_full = zv.key == "dr" ? z_int
        : Form("%.2f < %s < %.2f GeV (integrated)", azz->GetBinLowEdge(1), zv.sym.c_str(), azz->GetBinUpEdge(nz));

    if (zv.key == "dr") {
        // eps vs eta^pair in pT panels; eps vs pT in eta^pair panels (dR-integrated): the 2D
        // (pT, eta) cells, independent of the third axis -- drawn once, from the dR file.
        Figure F1 = base(outdir + "vs_pair_eta/pair_eta_in_pair_pt_bins_" + wp + "WP.png", "#eta^{pair}", false, z_int_full);
        for (int i = 1; i <= npt; ++i) {
            F1.panel_labels.push_back(PtLabel(apt, i));
            F1.curves.push_back({Eff(num, den, 'y', One(i), All(nz), Form("v1_%d", i))});
        }
        figs.push_back(F1);
        Figure F2 = base(outdir + "vs_pair_pt/pair_pt_in_pair_eta_bins_" + wp + "WP.png", "p_{T}^{pair} [GeV]", true, z_int_full);
        for (int j = 1; j <= neta; ++j) {
            F2.panel_labels.push_back(EtaLabel(aeta, j));
            F2.curves.push_back({Eff(num, den, 'x', One(j), All(nz), Form("v2_%d", j))});
        }
        figs.push_back(F2);
    }
    // eps vs Z in pT panels, and eps vs pT in Z panels: all |eta^pair| and per |eta^pair| group
    for (int g = -1; g < G.n; ++g) {
        const Sel eta_sel = (g < 0) ? All(neta) : G.ranges[g];
        const std::string tag = (g < 0) ? "abs_pair_eta_all" : AbsEtaTag(G, g);
        const std::string eta_txt = (g < 0) ? Form("|#eta^{pair}| < %.1f (integrated)", G.edges[G.n]) : AbsEtaLabel(G, g);
        Figure F3 = base(outdir + "vs_" + zv.key + "/" + zv.key + "_in_pair_pt_bins_" + tag + "_" + wp + "WP.png",
                         zv.sym + (zv.unit.empty() ? "" : " [GeV]"), false, eta_txt);
        for (int i = 1; i <= npt; ++i) {
            F3.panel_labels.push_back(PtLabel(apt, i));
            F3.curves.push_back({Eff(num, den, 'z', One(i), eta_sel, Form("v3_%d_%d", g, i))});
        }
        figs.push_back(F3);
        Figure F4 = base(outdir + "vs_pair_pt/pair_pt_in_" + zv.key + "_bins_" + tag + "_" + wp + "WP.png",
                         "p_{T}^{pair} [GeV]", true, eta_txt);
        for (int k = 1; k <= nz; ++k) {
            F4.panel_labels.push_back(ZLabel(zv, azz, k));
            F4.curves.push_back({Eff(num, den, 'x', eta_sel, One(k), Form("v4_%d_%d", g, k))});
        }
        figs.push_back(F4);
    }
    for (const auto& F : figs) {
        DrawFigure(F);
        const std::string xt = F.xtitle;
        const TAxis* xax = (xt.find("p_{T}") != std::string::npos) ? apt
                         : (xt.find("#eta") != std::string::npos) ? aeta : azz;
        WriteTable(tab, F, xax);
    }
    tab.close();
    std::cout << "[views] " << figs.size() << " figures -> " << outdir << std::endl;
}
