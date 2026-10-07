// =================================================================================================
// plot_pp24_fullsim_pair_bin_migration.cxx
//
// Bin-by-bin efficiency, stability and purity of the pp24 fullsim signal pairs (all truth OS, `_op`)
// on the canonical binnings, from the output of fill_pp24_fullsim_pair_bin_migration.C.
// docs/tracking/pair_bin_migration_bin_by_bin_feasibility.md (Physics Procedure §2, §3b).
//
//   efficiency  eps_i = N[T_i and R] / N[T_i]
//   stability   s_i   = N[T_i and R_i] / N[T_i and R]
//   purity      p_j   = N[T_j and R_j] / N[R_j]
//   outside     o_j   = N[R_j and not T] / N[R_j]
//   bin-by-bin  C_j   = N[T_j] / N[R_j]
//
// Output (under <sample>/plots/pp24_reco_effcy_plots/<wp>/signal_region/all_os/bin_migration/):
//   {pair_pt, pair_eta, minv, dr}/eff_stab_pur_<wp>WP.png   eps/s/p/o vs x + C_j
//   {pair_pt, pair_eta, minv, dr}/response_<wp>WP.png       P(reco bin | truth bin, selected)
//   cells_pt_eta_{minv,dr}/{purity,stability,efficiency}_<wp>WP.png   3D cells, one panel per pT bin
//   cells_compare/cell_purity_stability_<wp>WP.png          cumulative N[T] fraction vs cell p / s, m vs dR
//   bin_migration_values_<wp>WP.txt                         every plotted number
//
// Usage (from Analysis/plotting_codes/reco_effcy/bin_migration/):
//   root -l -b -q 'plot_pp24_fullsim_pair_bin_migration.cxx+(true)'    // Tight (nominal)
//   root -l -b -q 'plot_pp24_fullsim_pair_bin_migration.cxx+(false)'   // Medium (WP systematic)
// =================================================================================================

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "../../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../../RDFBasedHistFilling/CommonEffcyConfig.h"

namespace {

const double kNeffMin = 100.0;   // cells used in the summary min / fraction statistics
const double kNeffDraw = 10.0;   // 3D-cell panels: cells whose ratio DENOMINATOR has N_eff below this are not drawn
const std::vector<Color_t> kCols = {kBlack, kRed + 1, kBlue + 1, kGreen + 2};
const std::vector<Style_t> kMrks = {20, 21, 22, 24};

struct VarInfo {
    std::string key, title, unit, title_truth, title_reco;
    std::vector<double> edges;
    bool logx;
};

struct Axes {
    std::vector<double> pt, eta, minv, dr;
};

Axes CanonicalAxes()
{
    static const ParamsSet pms{};
    static const CommonEffcyConfig cfg{};
    Axes a;
    a.pt = pms.pair_pt_coarse_bins;
    a.eta.push_back(cfg.pair_eta_proj_ranges_coarse_incl_gap.front().first);
    for (const auto& r : cfg.pair_eta_proj_ranges_coarse_incl_gap) a.eta.push_back(r.second);
    a.minv = ParamsSet::SignalMinvRecoEffcyBins();
    a.dr.assign(cfg.dr_bins_edges_for_reco_effcy.begin(), cfg.dr_bins_edges_for_reco_effcy.end());
    return a;
}

template <class H>
H* Get(TFile& f, const std::string& n)
{
    auto* h = dynamic_cast<H*>(f.Get(n.c_str()));
    if (!h) throw std::runtime_error("plot_pp24_fullsim_pair_bin_migration: missing " + n);
    return h;
}

void CheckAxis(const TAxis* ax, const std::vector<double>& e, const std::string& what)
{
    if (ax->GetNbins() != (int)e.size() - 1)
        throw std::runtime_error("axis " + what + ": bin count differs from the canonical binning");
    for (int i = 1; i <= ax->GetNbins() + 1; ++i)
        if (std::fabs(ax->GetBinLowEdge(i) - e[i - 1]) > 1e-5)
            throw std::runtime_error("axis " + what + ": edge differs from the canonical binning");
}

// Ratio r = A/B with A a SUBSET of B (weighted binomial): var = [(1-2r) Sw2_A + r^2 Sw2_B] / B^2.
struct Val { double v = -1, e = 0, n = 0; bool ok = false; };   // n = N_eff of the denominator
Val SubsetRatio(double a, double ea2, double b, double eb2)
{
    Val r;
    if (b <= 0) return r;
    r.v = a / b;
    r.e = std::sqrt(std::max(0.0, ((1 - 2 * r.v) * ea2 + r.v * r.v * eb2) / (b * b)));
    r.n = eb2 > 0 ? b * b / eb2 : 0.0;
    r.ok = true;
    return r;
}
// 68 % Wilson score interval for a fraction r with N_eff = n (finite at r = 0 or 1, unlike the
// normal approximation): returns {lower error, upper error}.
std::pair<double, double> WilsonErr(double r, double n)
{
    const double z2 = 1.0, den = 1.0 + z2 / n;
    const double c = (r + z2 / (2 * n)) / den;
    const double h = std::sqrt(z2) / den * std::sqrt(r * (1 - r) / n + z2 / (4 * n * n));
    return {std::max(0.0, r - (c - h)), std::max(0.0, (c + h) - r)};
}
// Ratio of two overlapping populations treated as independent (an overestimate of the error).
Val PlainRatio(double a, double ea2, double b, double eb2)
{
    Val r;
    if (a <= 0 || b <= 0) return r;
    r.v = a / b;
    r.e = r.v * std::sqrt(ea2 / (a * a) + eb2 / (b * b));
    r.ok = true;
    return r;
}
double Neff(double s, double e2) { return e2 > 0 ? s * s / e2 : 0.0; }

// Per-bin quantities of one variable (1D) or of one 3D cell flattening.
struct Bins {
    int n = 0;
    std::vector<Val> eps, stab, pur, out, corr;
    std::vector<double> nT, neffT, neffTR, neffR;
};

Bins Compute(TFile& f, const std::string& key, const std::string& wp)
{
    auto* T  = Get<TH1D>(f, "h_" + key + "_T");
    auto* TR = Get<TH1D>(f, "h_" + key + "_TR_" + wp);
    auto* M  = Get<TH2D>(f, "h_" + key + "_M_" + wp);
    auto* R  = Get<TH1D>(f, "h_" + key + "_R_" + wp);
    auto* Ro = Get<TH1D>(f, "h_" + key + "_Rout_" + wp);
    Bins B;
    B.n = T->GetNbinsX();
    for (int i = 1; i <= B.n; ++i) {
        const double t = T->GetBinContent(i), t2 = std::pow(T->GetBinError(i), 2);
        const double tr = TR->GetBinContent(i), tr2 = std::pow(TR->GetBinError(i), 2);
        const double d = M->GetBinContent(i, i), d2 = std::pow(M->GetBinError(i, i), 2);
        const double r = R->GetBinContent(i), r2 = std::pow(R->GetBinError(i), 2);
        const double o = Ro->GetBinContent(i), o2 = std::pow(Ro->GetBinError(i), 2);
        B.eps.push_back(SubsetRatio(tr, tr2, t, t2));
        B.stab.push_back(SubsetRatio(d, d2, tr, tr2));
        B.pur.push_back(SubsetRatio(d, d2, r, r2));
        B.out.push_back(SubsetRatio(o, o2, r, r2));
        B.corr.push_back(PlainRatio(t, t2, r, r2));
        B.nT.push_back(t);
        B.neffT.push_back(Neff(t, t2));
        B.neffTR.push_back(Neff(tr, tr2));
        B.neffR.push_back(Neff(r, r2));
    }
    return B;
}

TH1D* ToHist(const std::vector<Val>& v, const std::vector<double>& edges, const std::string& name)
{
    auto* h = new TH1D(name.c_str(), "", (int)edges.size() - 1, edges.data());
    h->SetDirectory(nullptr);
    for (size_t i = 0; i < v.size(); ++i)
        if (v[i].ok) { h->SetBinContent(i + 1, v[i].v); h->SetBinError(i + 1, v[i].e); }
    return h;
}

std::string RangeLabel(const VarInfo& v, int i)
{
    const double lo = v.edges[i], hi = v.edges[i + 1];
    auto f = [&](double x) {
        char b[32];
        snprintf(b, sizeof b, v.key == "minv" ? "%.3f" : "%.1f", std::fabs(x));
        return std::string(x < -1e-9 ? "#minus" : "") + b;
    };
    return "[" + f(lo) + ", " + f(hi) + ")";
}
// "lo < title < hi unit" for panel titles
std::string RangeText(const VarInfo& v, int i, const std::string& unit)
{
    char b[96];
    snprintf(b, sizeof b, v.key == "minv" ? "%.3f < %s < %.3f%s" : "%.1f < %s < %.1f%s",
             v.edges[i], v.title.c_str(), v.edges[i + 1], unit.c_str());
    return b;
}

void Header(double x0, double y, double sz, const std::string& wp_text)
{
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextSize(sz);
    t.DrawLatex(x0, y, "Pythia 8 full simulation, #it{pp} 5.36 TeV, truth opposite-sign pairs");
    t.SetTextAlign(31);
    t.DrawLatex(0.97, y, wp_text.c_str());
}

// ---- 1D figure: eps, s, p, o vs x (top) and C_j (bottom) -------------------------------------
void Draw1D(const VarInfo& v, const Bins& B, const std::string& wp_text, const std::string& path)
{
    TCanvas c("c1d", "", 800, 760);
    TPad top("top", "", 0, 0.32, 1, 1), bot("bot", "", 0, 0, 1, 0.32);
    top.SetBottomMargin(0.02); top.SetTopMargin(0.17); top.SetLeftMargin(0.13); top.SetRightMargin(0.04);
    bot.SetTopMargin(0.03); bot.SetBottomMargin(0.34); bot.SetLeftMargin(0.13); bot.SetRightMargin(0.04);
    for (TPad* p : {&top, &bot}) { p->SetTicks(1, 1); p->SetLogx(v.logx); }
    top.Draw(); bot.Draw();
    top.cd();
    const std::vector<std::pair<const std::vector<Val>*, std::string>> series = {
        {&B.eps, "efficiency"}, {&B.stab, "stability"}, {&B.pur, "purity"},
        {&B.out, "outside-range fraction = N^{reco}(truth out of range) / N^{reco}"}};
    std::vector<TH1D*> hs;
    TLegend leg(0.13, 0.895, 0.96, 0.955), leg2(0.13, 0.835, 0.96, 0.895);   // row 1: eps/s/p; row 2: o
    leg.SetNColumns(3);
    for (TLegend* l : {&leg, &leg2}) { l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.042); }
    for (size_t s = 0; s < series.size(); ++s) {
        TH1D* h = ToHist(*series[s].first, v.edges, "h1d_" + std::to_string(s));
        h->SetLineColor(kCols[s]); h->SetMarkerColor(kCols[s]); h->SetMarkerStyle(kMrks[s]);
        h->SetMarkerSize(1.2); h->SetLineWidth(2); h->SetStats(0);
        h->GetYaxis()->SetRangeUser(0.0, 1.15);
        h->GetYaxis()->SetTitle("fraction");
        h->GetYaxis()->SetTitleSize(0.055); h->GetYaxis()->SetLabelSize(0.048); h->GetYaxis()->SetTitleOffset(1.0);
        h->GetXaxis()->SetLabelSize(0); h->GetXaxis()->SetTitleSize(0);
        h->Draw(s == 0 ? "PE" : "PE SAME");
        (s < 3 ? leg : leg2).AddEntry(h, series[s].second.c_str(), "lp");
        hs.push_back(h);
    }
    leg.Draw();
    leg2.Draw();
    Header(0.13, 0.965, 0.042, wp_text);
    bot.cd();
    TH1D* hc = ToHist(B.corr, v.edges, "h1d_corr");
    double cmax = 0, cmin = 1e9;
    for (const auto& x : B.corr) if (x.ok) { cmax = std::max(cmax, x.v + x.e); cmin = std::min(cmin, x.v - x.e); }
    hc->SetStats(0); hc->SetMarkerStyle(20); hc->SetLineWidth(2);
    hc->GetYaxis()->SetRangeUser(std::floor(cmin * 5) / 5 - 0.1, std::ceil(cmax * 5) / 5 + 0.1);
    hc->GetYaxis()->SetNdivisions(505);
    hc->GetYaxis()->SetTitle("C_{j} = N_{j}^{truth} / N_{j}^{reco}");
    hc->GetYaxis()->SetTitleSize(0.10); hc->GetYaxis()->SetLabelSize(0.10); hc->GetYaxis()->SetTitleOffset(0.55);
    hc->GetXaxis()->SetTitle((v.title + v.unit).c_str());
    hc->GetXaxis()->SetTitleSize(0.12); hc->GetXaxis()->SetLabelSize(0.10); hc->GetXaxis()->SetTitleOffset(1.1);
    hc->GetXaxis()->SetNoExponent(); hc->GetXaxis()->SetMoreLogLabels(v.logx);
    hc->Draw("PE");
    gSystem->mkdir(gSystem->DirName(path.c_str()), kTRUE);
    c.SaveAs(path.c_str());
    for (auto* h : hs) delete h;
    delete hc;
}

// ---- response: P(reco bin j | truth bin i, selected) -------------------------------------------
void DrawResponse(TFile& f, const VarInfo& v, const std::string& wp, const std::string& wp_text,
                  const std::string& path)
{
    auto* M  = Get<TH2D>(f, "h_" + v.key + "_M_" + wp);
    auto* TR = Get<TH1D>(f, "h_" + v.key + "_TR_" + wp);
    const int n = M->GetNbinsX();
    TH2D P("P", "", n, 0, n, n, 0, n);
    P.SetDirectory(nullptr);
    for (int i = 1; i <= n; ++i) {
        P.GetXaxis()->SetBinLabel(i, RangeLabel(v, i - 1).c_str());
        P.GetYaxis()->SetBinLabel(i, RangeLabel(v, i - 1).c_str());
        const double tr = TR->GetBinContent(i);
        for (int j = 1; j <= n; ++j)
            if (tr > 0) P.SetBinContent(i, j, M->GetBinContent(i, j) / tr);
    }
    TCanvas c("cresp", "", 900, 820);
    c.SetLeftMargin(0.19); c.SetBottomMargin(0.17); c.SetRightMargin(0.14); c.SetTopMargin(0.08);
    gStyle->SetPaintTextFormat(".3g");
    P.SetStats(0);
    P.SetMinimum(0.0); P.SetMaximum(1.0);
    P.SetMarkerSize(n > 8 ? 0.85 : 1.1);
    P.GetXaxis()->SetTitle((v.title_truth + v.unit).c_str());
    P.GetYaxis()->SetTitle((v.title_reco + v.unit).c_str());
    P.GetXaxis()->SetTitleOffset(1.75); P.GetYaxis()->SetTitleOffset(2.3);
    P.GetXaxis()->SetLabelSize(0.028); P.GetYaxis()->SetLabelSize(0.028);
    P.GetZaxis()->SetTitle("P(reco bin | truth bin, selected)");
    P.GetZaxis()->SetTitleOffset(1.3);
    P.Draw("COLZ TEXT");
    Header(0.19, 0.94, 0.028, wp_text);
    gSystem->mkdir(gSystem->DirName(path.c_str()), kTRUE);
    c.SaveAs(path.c_str());
}

// ---- 3D cells: one panel per pT bin, x = eta^pair, series = Z bins ----------------------------
void DrawCells(const Bins& B, const std::vector<Val> Bins::*q, const std::string& ytitle,
               const VarInfo& pt, const VarInfo& eta, const VarInfo& z, const std::string& wp_text,
               const std::string& path)
{
    const int npt = pt.edges.size() - 1, neta = eta.edges.size() - 1, nz = z.edges.size() - 1;
    // y range common to all panels: 0.1 below the smallest drawn value, rounded down, at most 0.8
    double vmin = 1.0;
    for (const auto& x : B.*q) if (x.ok && x.n >= kNeffDraw) vmin = std::min(vmin, x.v);
    const double ylo = std::max(0.0, std::min(0.8, std::floor(vmin * 10) / 10 - 0.1));
    const int ncol = (int)std::ceil(std::sqrt((double)npt)), nrow = (int)std::ceil((double)npt / ncol);
    const int hstrip = 110;
    TCanvas c("ccells", "", 430 * ncol, 340 * nrow + hstrip);
    const double fh = (double)hstrip / (340.0 * nrow + hstrip);
    TPad* head = new TPad("head", "", 0, 1 - fh, 1, 1);
    TPad* body = new TPad("body", "", 0, 0, 1, 1 - fh);
    head->Draw(); body->Draw();
    std::vector<TGraphAsymmErrors*> gs;
    std::vector<TH1D*> frames;
    body->cd();
    body->Divide(ncol, nrow);
    for (int ip = 0; ip < npt; ++ip) {
        body->cd(ip + 1);
        gPad->SetTicks(1, 1); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.12);
        auto* fr = new TH1D(Form("fr%d", ip), "", neta, eta.edges.data());
        fr->SetDirectory(nullptr); fr->SetStats(0);
        fr->GetYaxis()->SetRangeUser(ylo, 1.08);
        fr->GetXaxis()->SetTitle(eta.title.c_str());
        fr->GetYaxis()->SetTitle(ytitle.c_str());
        fr->GetXaxis()->SetTitleSize(0.052); fr->GetXaxis()->SetLabelSize(0.044);
        fr->GetYaxis()->SetTitleSize(0.052); fr->GetYaxis()->SetLabelSize(0.044); fr->GetYaxis()->SetTitleOffset(1.35);
        fr->Draw("AXIS");
        frames.push_back(fr);
        for (int iz = 0; iz < nz; ++iz) {
            auto* g = new TGraphAsymmErrors();
            for (int ie = 0; ie < neta; ++ie) {
                const int k = (ip * neta + ie) * nz + iz;
                const Val& x = (B.*q)[k];
                if (!x.ok || x.n < kNeffDraw) continue;
                const double w = eta.edges[ie + 1] - eta.edges[ie];
                const double xc = 0.5 * (eta.edges[ie] + eta.edges[ie + 1]) + (iz - 0.5 * (nz - 1)) * 0.16 * w;
                const int n = g->GetN();
                g->SetPoint(n, xc, x.v);
                const auto [elo, ehi] = WilsonErr(x.v, x.n);
                g->SetPointError(n, 0, 0, elo, ehi);
            }
            g->SetMarkerColor(kCols[iz % kCols.size()]); g->SetLineColor(kCols[iz % kCols.size()]);
            g->SetMarkerStyle(kMrks[iz % kMrks.size()]); g->SetMarkerSize(0.9); g->SetLineWidth(1);
            if (g->GetN() > 0) g->Draw("P SAME");
            gs.push_back(g);
        }
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.052);
        t.DrawLatex(0.17, 0.915, RangeText(pt, ip, " GeV").c_str());
    }
    head->cd();
    Header(0.02, 0.80, 0.24 * 110.0 / hstrip * (ncol <= 2 ? 0.85 : 1.0), wp_text);
    TLegend* leg = new TLegend(0.02, 0.05, 0.98, 0.55);
    leg->SetNColumns(nz); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.22);
    for (int iz = 0; iz < nz; ++iz)
        leg->AddEntry(gs[iz], RangeText(z, iz, z.unit).c_str(), "p");
    leg->Draw();
    gSystem->mkdir(gSystem->DirName(path.c_str()), kTRUE);
    c.SaveAs(path.c_str());
    for (auto* g : gs) delete g;
    for (auto* h : frames) delete h;
}

// ---- cumulative N[T]-weighted fraction of cells with purity / stability below x, m vs dR -------
// F(x) = sum_{cells: q < x} N[T] / sum_{cells: q defined} N[T]. Cumulative, so nothing is hidden by
// the axis range: F at the left edge is the fraction below it. LOG y (documented exception to the
// linear default): the m-vs-dR decision rests on the low-purity TAILS, which are 1e-4..1e-2 of
// N[T] next to F -> 1 at x = 1 -- unreadable on a linear axis.
void DrawCellCompare(const Bins& Bm, const Bins& Bd, const std::string& wp_text, const std::string& path)
{
    TCanvas c("ccmp", "", 1300, 600);
    TPad* head = new TPad("head", "", 0, 0.86, 1, 1);
    TPad* body = new TPad("body", "", 0, 0, 1, 0.86);
    head->Draw(); body->Draw();
    body->cd();
    body->Divide(2, 1);
    std::vector<TH1D*> keep;
    const std::vector<std::pair<const std::vector<Val> Bins::*, std::string>> qs = {
        {&Bins::pur, "purity"}, {&Bins::stab, "stability"}};
    const double xlo = 0.5, xhi = 1.0;
    const int nb = 100;
    TH1D* leg_m = nullptr; TH1D* leg_d = nullptr;
    for (int iq = 0; iq < 2; ++iq) {
        body->cd(iq + 1);
        gPad->SetTicks(1, 1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
        gPad->SetLogy(true);   // see the function comment
        auto mk = [&](const Bins& B, const std::string& n, Color_t col, Style_t ls) {
            auto* h = new TH1D(n.c_str(), "", nb, xlo, xhi);
            h->SetDirectory(nullptr);
            double tot = 0;
            for (int k = 0; k < B.n; ++k)
                if ((B.*(qs[iq].first))[k].ok) tot += B.nT[k];
            for (int b = 1; b <= nb; ++b) {
                const double x = h->GetXaxis()->GetBinUpEdge(b);
                double below = 0;
                for (int k = 0; k < B.n; ++k) {
                    const Val& v = (B.*(qs[iq].first))[k];
                    if (v.ok && v.v < x) below += B.nT[k];
                }
                h->SetBinContent(b, below / tot);
            }
            h->SetLineColor(col); h->SetLineWidth(2); h->SetLineStyle(ls); h->SetStats(0);
            keep.push_back(h);
            return h;
        };
        TH1D* hm = mk(Bm, Form("hm%d", iq), kRed + 1, 1);
        TH1D* hd = mk(Bd, Form("hd%d", iq), kBlue + 1, 7);
        // floor below the smallest nonzero F of either curve, so no part of either curve leaves the frame
        double fmin = 1.0;
        for (TH1D* h : {hm, hd})
            for (int b = 1; b <= nb; ++b)
                if (h->GetBinContent(b) > 0) fmin = std::min(fmin, h->GetBinContent(b));
        hm->GetYaxis()->SetRangeUser(std::pow(10.0, std::floor(std::log10(fmin))), 2.0);
        hm->GetXaxis()->SetTitle(("3D-cell " + qs[iq].second + " x").c_str());
        hm->GetYaxis()->SetTitle(("fraction of N^{truth} in cells with " + qs[iq].second + " < x").c_str());
        hm->GetXaxis()->SetTitleSize(0.05); hm->GetYaxis()->SetTitleSize(0.045);
        hm->GetXaxis()->SetLabelSize(0.045); hm->GetYaxis()->SetLabelSize(0.045); hm->GetYaxis()->SetTitleOffset(1.45);
        hm->Draw("HIST L");
        hd->Draw("HIST L SAME");
        leg_m = hm; leg_d = hd;
    }
    head->cd();
    Header(0.02, 0.62, 0.30, wp_text);
    TLegend* leg = new TLegend(0.02, 0.0, 0.98, 0.42);
    leg->SetNColumns(2); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.28);
    leg->AddEntry(leg_m, "(p_{T}^{pair}, #eta^{pair}, m_{#mu#mu}) cells", "l");
    leg->AddEntry(leg_d, "(p_{T}^{pair}, #eta^{pair}, #DeltaR) cells", "l");
    leg->Draw();
    gSystem->mkdir(gSystem->DirName(path.c_str()), kTRUE);
    c.SaveAs(path.c_str());
    for (auto* h : keep) delete h;
}

void WriteBins(std::ofstream& o, const std::string& name, const Bins& B,
               const std::vector<std::string>& labels)
{
    o << "\n## " << name << "\n";
    o << "bin | N_eff[T] | eps | stability | purity | outside | C\n";
    for (int i = 0; i < B.n; ++i) {
        auto s = [](const Val& v) { char b[48]; if (!v.ok) return std::string("n/a");
                                    snprintf(b, sizeof b, "%.4f +- %.4f", v.v, v.e); return std::string(b); };
        o << labels[i] << " | " << Form("%.0f", B.neffT[i]) << " | " << s(B.eps[i]) << " | " << s(B.stab[i])
          << " | " << s(B.pur[i]) << " | " << s(B.out[i]) << " | " << s(B.corr[i]) << "\n";
    }
}

void WriteCellSummary(std::ofstream& o, const std::string& name, const Bins& B)
{
    double sT = 0, sTs = 0, sTp = 0, sTstab = 0, sTpur = 0, wlow_p = 0, wlow_s = 0, wlow_p5 = 0;
    double minp = 2, mins = 2; int nused = 0, nall = 0;
    for (int k = 0; k < B.n; ++k) {
        if (B.nT[k] <= 0) continue;
        ++nall;
        sT += B.nT[k];
        // each mean over the cells where that quantity is defined
        if (B.stab[k].ok) { sTstab += B.nT[k] * B.stab[k].v; sTs += B.nT[k]; }
        if (B.pur[k].ok)  { sTpur  += B.nT[k] * B.pur[k].v;  sTp += B.nT[k]; }
        if (B.pur[k].ok && B.pur[k].v < 0.8)  wlow_p += B.nT[k];
        if (B.pur[k].ok && B.pur[k].v < 0.5)  wlow_p5 += B.nT[k];
        if (B.stab[k].ok && B.stab[k].v < 0.8) wlow_s += B.nT[k];
        if (B.neffTR[k] >= kNeffMin && B.neffR[k] >= kNeffMin) {
            ++nused;
            minp = std::min(minp, B.pur[k].v);
            mins = std::min(mins, B.stab[k].v);
        }
    }
    o << "\n## " << name << " -- 3D-cell summary\n"
      << "cells with N[T] > 0: " << nall << "\n"
      << Form("N[T]-weighted mean stability %.4f, purity %.4f\n", sTstab / sTs, sTpur / sTp)
      << Form("cells with N_eff >= %.0f (TR and R): %d; min stability %.4f, min purity %.4f\n", kNeffMin, nused, mins, minp)
      << Form("fraction of N[T] in cells with purity < 0.8: %.4f, < 0.5: %.4f; stability < 0.8: %.4f\n",
              wlow_p / sT, wlow_p5 / sT, wlow_s / sT);
}

} // namespace

void plot_pp24_fullsim_pair_bin_migration(bool tight_WP = true)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    const std::string wp = tight_WP ? "tight" : "medium";
    const std::string wp_text = tight_WP ? "Tight WP" : "Medium WP";

    const std::string sample_dir = FullSimSampleInputDir(FullSimSampleType::pp, /*is_test_sample=*/false);
    const std::string in_file = FullSimRecoEffDir(sample_dir) + "pair_bin_migration_pp24_full.root";
    const std::string out = FullSimPlotsDir(sample_dir) + "pp24_reco_effcy_plots/" + wp
                          + "/signal_region/all_os/bin_migration/";
    TFile f(in_file.c_str(), "READ");
    if (f.IsZombie()) throw std::runtime_error("cannot open " + in_file);

    const Axes A = CanonicalAxes();
    const VarInfo vpt  {"pair_pt",  "p_{T}^{pair}", " [GeV]", "p_{T}^{pair,truth}", "p_{T}^{pair,reco}", A.pt,   true};
    const VarInfo veta {"pair_eta", "#eta^{pair}",  "",       "#eta^{pair,truth}",  "#eta^{pair,reco}",  A.eta,  false};
    const VarInfo vm   {"minv",     "m_{#mu#mu}",   " [GeV]", "m_{#mu#mu}^{truth}", "m_{#mu#mu}^{reco}", A.minv, false};
    const VarInfo vdr  {"dr",       "#DeltaR",      "",       "#DeltaR^{truth}",    "#DeltaR^{reco}",    A.dr,   false};
    VarInfo vm_lbl = vm;
    vm_lbl.unit = " GeV";   // legend form of the mass label

    std::ofstream tab;
    gSystem->mkdir(out.c_str(), kTRUE);
    tab.open(out + "bin_migration_values_" + wp + "WP.txt");
    tab << "# pp24 fullsim, all truth OS pairs, " << wp_text << "; input " << in_file << "\n"
        << "# eps=N[T_i&R]/N[T_i]  stability=N[T_i&R_i]/N[T_i&R]  purity=N[T_j&R_j]/N[R_j]  "
           "outside=N[R_j&!T]/N[R_j]  C=N[T_j]/N[R_j] (C error: T,R treated as independent)\n";

    for (const VarInfo* v : {&vpt, &veta, &vm, &vdr}) {
        CheckAxis(Get<TH1D>(f, "h_" + v->key + "_T")->GetXaxis(), v->edges, v->key);
        const Bins B = Compute(f, v->key, wp);
        Draw1D(*v, B, wp_text, out + v->key + "/eff_stab_pur_" + wp + "WP.png");
        DrawResponse(f, *v, wp, wp_text, out + v->key + "/response_" + wp + "WP.png");
        std::vector<std::string> lab;
        for (int i = 0; i < B.n; ++i) lab.push_back(Form("[%.3f, %.3f)", v->edges[i], v->edges[i + 1]));
        WriteBins(tab, v->key, B, lab);
    }

    const Bins Bm = Compute(f, "cell_minv", wp);
    const Bins Bd = Compute(f, "cell_dr", wp);
    for (const auto& [B, z, tag] : {std::make_tuple(&Bm, (const VarInfo*)&vm_lbl, std::string("minv")),
                                    std::make_tuple(&Bd, &vdr, std::string("dr"))}) {
        if (B->n != (int)((vpt.edges.size() - 1) * (veta.edges.size() - 1) * (z->edges.size() - 1)))
            throw std::runtime_error("cells_pt_eta_" + tag + ": cell count differs from the canonical axes (stale input?)");
        const std::string d = out + "cells_pt_eta_" + tag + "/";
        DrawCells(*B, &Bins::pur,  "purity",            vpt, veta, *z, wp_text, d + "purity_" + wp + "WP.png");
        DrawCells(*B, &Bins::stab, "stability",         vpt, veta, *z, wp_text, d + "stability_" + wp + "WP.png");
        DrawCells(*B, &Bins::eps,  "efficiency", vpt, veta, *z, wp_text, d + "efficiency_" + wp + "WP.png");
        const int neta = veta.edges.size() - 1, nz = z->edges.size() - 1;
        std::vector<std::string> lab;
        for (int k = 0; k < B->n; ++k) {
            const int ip = k / (neta * nz), ie = (k / nz) % neta, iz = k % nz;
            lab.push_back(Form("pT[%.1f,%.1f) eta[%.1f,%.1f) %s[%.3f,%.3f)", vpt.edges[ip], vpt.edges[ip + 1],
                               veta.edges[ie], veta.edges[ie + 1], tag.c_str(), z->edges[iz], z->edges[iz + 1]));
        }
        WriteCellSummary(tab, "cells_pt_eta_" + tag, *B);
        WriteBins(tab, "cells_pt_eta_" + tag, *B, lab);
    }
    DrawCellCompare(Bm, Bd, wp_text, out + "cells_compare/cell_purity_stability_" + wp + "WP.png");
    tab.close();
    std::cout << "[done] " << out << "\n";
}
