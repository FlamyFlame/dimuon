// =================================================================================================
// plot_pp24_fullsim_pair_reco_eff_reco_binned.cxx
//
// Views of the RECO-BINNED pair reco efficiency eps^reco = N[R](reco) / N[T](truth) on the canonical
// (pair pT, eta^pair, m_mumu) cells, all truth OS pairs, with the TRUTH-binned eps = N[T and R](truth)
// / N[T](truth) on the same axes overlaid for comparison
// (docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md §3b).
//
// Every point is a ratio of SUMMED counts over its set of cells (reco_binned_cells.h) -- never a
// projection of a ratio. The single-cell values are cross-checked against the applied product
// (pair_reco_eff_reco_binned_pp24_full.root); the macro THROWS on a mismatch.
//
// Figures (dir .../plots/pp24_reco_effcy_plots/<wp>/signal_region/all_os/reco_binned/):
//   vs_pair_eta/pair_eta_in_pair_pt_bins_<wp>WP.png   (m-integrated)
//   vs_pair_pt/pair_pt_in_pair_eta_bins_<wp>WP.png    (m-integrated)
//   vs_minv/minv_in_pair_pt_bins_<wp>WP.png           (eta^pair-integrated)
//   reco_binned_values_<wp>WP.txt
//
// Usage (from Analysis/plotting_codes/reco_effcy/reco_binned/):
//   root -l -b -q 'plot_pp24_fullsim_pair_reco_eff_reco_binned.cxx+(true)'    // Tight (nominal)
//   root -l -b -q 'plot_pp24_fullsim_pair_reco_eff_reco_binned.cxx+(false)'   // Medium (WP systematic)
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
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include "../../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../../Utilities/PairRecoEffRecoBinnedEvaluator.h"
#include "reco_binned_cells.h"

namespace {

struct Panel {
    std::string label;
    std::vector<CellRatio> pts;   // one per x bin
};
struct Fig {
    std::string path, xtitle, extra;
    std::vector<double> xedges;
    bool logx = false;
    std::vector<Panel> panels;
};

std::string PtLabel(const std::vector<double>& e, int i)
{ return Form("%.4g < p_{T}^{pair} < %.4g GeV", e[i], e[i + 1]); }
std::string EtaLabel(const std::vector<double>& e, int j)
{ return Form("%.1f < #eta^{pair} < %.1f", e[j], e[j + 1]); }
std::string Ascii(std::string s)
{
    const std::vector<std::pair<std::string, std::string>> rep = {
        {"#eta^{pair}", "eta_pair"}, {"p_{T}^{pair}", "pT_pair"}, {"m_{#mu#mu}", "m_mumu"}};
    for (const auto& r : rep)
        for (size_t p; (p = s.find(r.first)) != std::string::npos;) s.replace(p, r.first.size(), r.second);
    return s;
}

// x position of series s (0 reco-binned, 1 truth-binned) inside bin b: offset +-12 % of the bin
// width (multiplicative on log x) so the two error bars do not hide each other.
double XPos(const Fig& F, int b, int s)
{
    const double lo = F.xedges[b], hi = F.xedges[b + 1], f = (s == 0) ? 0.38 : 0.62;
    return F.logx ? std::exp(std::log(lo) + f * (std::log(hi) - std::log(lo))) : lo + f * (hi - lo);
}

void Draw(const Fig& F, const std::string& header, const std::string& wp_text)
{
    const int np = F.panels.size();
    const int ncol = (np <= 3) ? np : (int)std::ceil(std::sqrt((double)np));   // nrows >= ncols
    const int nrow = (int)std::ceil((double)np / ncol);
    const int hstrip = 150;
    TCanvas c("c_fig", "", 430 * ncol, 340 * nrow + hstrip);
    const double fh = (double)hstrip / (340.0 * nrow + hstrip);
    TPad* head = new TPad("head", "", 0.0, 1.0 - fh, 1.0, 1.0);
    TPad* body = new TPad("body", "", 0.0, 0.0, 1.0, 1.0 - fh);
    head->SetFillStyle(0); body->SetFillStyle(0);
    head->Draw(); body->Draw();

    // graphs first (legend needs them)
    std::vector<std::vector<TGraphErrors*>> g(np, std::vector<TGraphErrors*>(2, nullptr));
    double ymax = 1.1;
    for (int p = 0; p < np; ++p)
        for (int s = 0; s < 2; ++s) {
            auto* gr = new TGraphErrors();
            for (size_t b = 0; b < F.panels[p].pts.size(); ++b) {
                const CellRatio& r = F.panels[p].pts[b];
                const bool ok = (s == 0) ? r.measured : r.truth_measured;
                if (!ok) continue;
                const double v = (s == 0) ? r.eps : r.eps_truth, e = (s == 0) ? r.err : r.err_truth;
                const int n = gr->GetN();
                gr->SetPoint(n, XPos(F, b, s), v);
                gr->SetPointError(n, 0.0, e);
                ymax = std::max(ymax, 1.08 * (v + e));
            }
            gr->SetMarkerStyle(s == 0 ? 20 : 24);
            gr->SetMarkerColor(s == 0 ? kBlack : kRed + 1);
            gr->SetLineColor(s == 0 ? kBlack : kRed + 1);
            gr->SetLineWidth(2);
            gr->SetMarkerSize(1.1);
            g[p][s] = gr;
        }

    head->cd();
    {
        const double sz = (ncol <= 2 ? 0.20 : 0.24) * 110.0 / hstrip;
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(sz);
        double y = 0.80;
        t.DrawLatex(0.02, y, header.c_str());
        t.SetTextAlign(31); t.DrawLatex(0.98, y, wp_text.c_str()); t.SetTextAlign(11);
        y -= 1.35 * sz;
        t.DrawLatex(0.02, y, "All opposite-sign pairs, truth signal region; N_{reco} incl. detector-gap cuts");
        y -= 1.35 * sz;
        t.DrawLatex(0.02, y, F.extra.c_str());
        TLegend* leg = new TLegend(0.02, 0.02, 0.98, y - 0.6 * sz);
        leg->SetNColumns(2); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(sz);
        leg->AddEntry(g[0][0], "reco-binned: N_{reco}(reco bin) / N_{gen}(truth bin)", "lp");
        leg->AddEntry(g[0][1], "truth-binned: N_{reco}(truth bin) / N_{gen}(truth bin)", "lp");
        leg->Draw();
    }
    body->cd();
    body->Divide(ncol, nrow);
    for (int p = 0; p < np; ++p) {
        body->cd(p + 1);
        gPad->SetTicks(1, 1); gPad->SetLogx(F.logx);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.12);
        auto* fr = new TH1D(Form("fr%d", p), "", (int)F.xedges.size() - 1, F.xedges.data());
        fr->SetStats(0);
        fr->GetXaxis()->SetTitle(F.xtitle.c_str());
        fr->GetXaxis()->SetNoExponent(); fr->GetXaxis()->SetMoreLogLabels(F.logx);
        fr->GetYaxis()->SetTitle("#varepsilon_{reco}^{pair}");
        fr->GetYaxis()->SetRangeUser(0.0, ymax);
        fr->GetXaxis()->SetTitleSize(0.052); fr->GetXaxis()->SetLabelSize(0.044);
        fr->GetYaxis()->SetTitleSize(0.052); fr->GetYaxis()->SetLabelSize(0.044); fr->GetYaxis()->SetTitleOffset(1.2);
        fr->Draw("AXIS");
        TLine* one = new TLine(F.xedges.front(), 1.0, F.xedges.back(), 1.0);
        one->SetLineStyle(2); one->SetLineColor(kGray + 1); one->Draw();
        g[p][1]->Draw("P SAME");
        g[p][0]->Draw("P SAME");
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.052);
        t.DrawLatex(0.17, 0.915, F.panels[p].label.c_str());
    }
    gSystem->mkdir(gSystem->DirName(F.path.c_str()), kTRUE);
    c.SaveAs(F.path.c_str());
}

void Table(std::ofstream& tab, const Fig& F)
{
    tab << "\n### " << gSystem->BaseName(F.path.c_str()) << "  (x = " << Ascii(F.xtitle) << "; " << Ascii(F.extra) << ")\n";
    for (const auto& P : F.panels) {
        tab << "  [" << Ascii(P.label) << "]\n";
        for (size_t b = 0; b < P.pts.size(); ++b) {
            const CellRatio& r = P.pts[b];
            const std::string sr = r.measured ? std::string(Form("%.4f+-%.4f", r.eps, r.err)) : std::string("-");
            const std::string st = r.truth_measured ? std::string(Form("%.4f+-%.4f", r.eps_truth, r.err_truth)) : std::string("-");
            tab << Form("    [%g,%g]  reco-binned ", F.xedges[b], F.xedges[b + 1]) << sr << "   truth-binned " << st
                << Form("   N[R]=%.4g N[T]=%.4g A=%.4g Neff[T]=%.0f\n", r.nR, r.nT, r.nA, r.NeffT());
        }
    }
}

}  // namespace

void plot_pp24_fullsim_pair_reco_eff_reco_binned(bool tight_WP = true)
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    using EV = PairRecoEffRecoBinnedEvaluator;
    const std::vector<double> ept = EV::PtEdges(), eeta = EV::EtaEdges(), em = EV::MinvEdges();
    const int npt = ept.size() - 1, neta = eeta.size() - 1, nm = em.size() - 1;
    const std::string wp = tight_WP ? "tight" : "medium";
    const std::string wp_text = tight_WP ? "Tight muons" : "Medium muons";
    const std::string header = "Pythia8 full simulation, pp 2024 conditions";

    const std::string sample_dir = FullSimSampleInputDir(FullSimSampleType::pp, false);
    const std::string mig_path = FullSimRecoEffDir(sample_dir) + "pair_bin_migration_pp24_full.root";
    const std::string prod_path = FullSimRecoEffDir(sample_dir) + "pair_reco_eff_reco_binned_pp24_full.root";
    const std::string outdir = FullSimPlotsDir(sample_dir) + "pp24_reco_effcy_plots/" + wp
                             + "/signal_region/all_os/reco_binned/";

    std::unique_ptr<TFile> fm(TFile::Open(mig_path.c_str(), "READ"));
    if (!fm || fm->IsZombie()) throw std::runtime_error("cannot open " + mig_path);
    const MigrationCells C = LoadMigrationCells(fm.get(), wp, npt, neta, nm);

    // Cross-check against the applied product (marker, axes, every cell).
    PairRecoEffRecoBinnedEvaluator ev;
    ev.Load(prod_path, tight_WP);
    double maxdiff = 0;
    for (int i = 0; i < npt; ++i) for (int j = 0; j < neta; ++j) {
        for (int k = 0; k < nm; ++k) {
            const CellRatio r = RatioOverCells(C, {C.Cell(i, j, k)});
            maxdiff = std::max(maxdiff, std::fabs((r.measured ? r.eps : 0.) - ev.h_eff->GetBinContent(i + 1, j + 1, k + 1)));
        }
        const CellRatio r2 = RatioOverCells(C, CellsIn(C, i, i, j, j, 0, nm - 1));
        maxdiff = std::max(maxdiff, std::fabs((r2.measured ? r2.eps : 0.) - ev.h_eff2d->GetBinContent(i + 1, j + 1)));
    }
    if (maxdiff > 1e-9) throw std::runtime_error("plotter cells differ from the product by " + std::to_string(maxdiff));
    std::cout << "[check] plotted cells == product (max |diff| " << maxdiff << ")\n";

    std::vector<Fig> figs;
    {
        Fig F; F.path = outdir + "vs_pair_eta/pair_eta_in_pair_pt_bins_" + wp + "WP.png";
        F.xtitle = "#eta^{pair}"; F.xedges = eeta;
        F.extra = Form("%g < m_{#mu#mu} < %g GeV (integrated)", em.front(), em.back());
        for (int i = 0; i < npt; ++i) {
            Panel P; P.label = PtLabel(ept, i);
            for (int j = 0; j < neta; ++j) P.pts.push_back(RatioOverCells(C, CellsIn(C, i, i, j, j, 0, nm - 1)));
            F.panels.push_back(P);
        }
        figs.push_back(F);
    }
    {
        Fig F; F.path = outdir + "vs_pair_pt/pair_pt_in_pair_eta_bins_" + wp + "WP.png";
        F.xtitle = "p_{T}^{pair} [GeV]"; F.xedges = ept; F.logx = true;
        F.extra = Form("%g < m_{#mu#mu} < %g GeV (integrated)", em.front(), em.back());
        for (int j = 0; j < neta; ++j) {
            Panel P; P.label = EtaLabel(eeta, j);
            for (int i = 0; i < npt; ++i) P.pts.push_back(RatioOverCells(C, CellsIn(C, i, i, j, j, 0, nm - 1)));
            F.panels.push_back(P);
        }
        figs.push_back(F);
    }
    {
        Fig F; F.path = outdir + "vs_minv/minv_in_pair_pt_bins_" + wp + "WP.png";
        F.xtitle = "m_{#mu#mu} [GeV]"; F.xedges = em;
        F.extra = Form("%.1f < #eta^{pair} < %.1f (integrated)", eeta.front(), eeta.back());
        for (int i = 0; i < npt; ++i) {
            Panel P; P.label = PtLabel(ept, i);
            for (int k = 0; k < nm; ++k) P.pts.push_back(RatioOverCells(C, CellsIn(C, i, i, 0, neta - 1, k, k)));
            F.panels.push_back(P);
        }
        figs.push_back(F);
    }

    gSystem->mkdir(outdir.c_str(), kTRUE);
    std::ofstream tab(outdir + "reco_binned_values_" + wp + "WP.txt");
    std::vector<int> all;
    for (int c = 0; c < npt * neta * nm; ++c) all.push_back(c);
    const CellRatio ri = RatioOverCells(C, all);
    int n_gt1 = 0;
    for (int c : all) { const CellRatio r = RatioOverCells(C, {c}); if (r.measured && r.eps > 1) ++n_gt1; }
    tab << "Reco-binned pair reco efficiency eps^reco_S = N[R_S](reco) / N[T_S](truth) (pp24 fullsim FULL, " << wp_text
        << ", all truth OS pairs), ratio of summed counts over each cell set S; truth-binned eps = N[T_S and R](truth)/N[T_S].\n"
        << "Input " << mig_path << "\nProduct " << prod_path << "\n"
        << Form("Inclusive: reco-binned %.5f+-%.5f, truth-binned %.5f+-%.5f (N[R]=%.5g, N[T]=%.5g, A=%.5g)\n",
                ri.eps, ri.err, ri.eps_truth, ri.err_truth, ri.nR, ri.nT, ri.nA)
        << "3D cells with eps^reco > 1 (not capped): " << n_gt1 << "\n";
    for (const auto& F : figs) { Draw(F, header, wp_text); Table(tab, F); }
    tab.close();
    std::cout << "[done] " << figs.size() << " figures -> " << outdir << std::endl;
}
