// =================================================================================================
// plot_ss_fraction_vs_pair_pt.cxx
//
// DATA raw same-sign (SS) and opposite-sign (OS) pair counts in the single-b signal region vs pair
// pT, and the SS fraction f_SS = N_SS / (N_SS + N_OS). 2x2 canvas: pp24 (left) | Pb+Pb 50-80 %,
// 2023+2024+2025+2026 summed (right); raw counts (top) | f_SS (bottom).
// docs/tracking/ss_fraction_vs_pair_pt_data.md.
//
//   input : crossx hist-filling outputs (RDFBasedHistFillingPP/PbPb.cxx), UNWEIGHTED counts from
//           the same RDF nodes as the SS/OS crossx yields (Tight WP + signal_cuts, i.e. the signal
//           mass window):
//     pp   : h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts         (OS)
//            h2d_ss_counts_pair_pt_pair_eta_binned_w_signal_cuts      (SS)
//     PbPb : h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_ctr50_80_counts   (OS)
//            h2d_ss_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_ctr50_80_counts   (SS)
//   Pair-pT axis = the histograms' own (ParamsSet::pT_bins_150, the default crossx axis);
//   projected over all pair-eta bins (the signal region already requires |eta^pair| < 2.2).
//   f_SS error = Clopper-Pearson 68.3 % interval.
//
// Output: plots/single_b_analysis/ss_fraction/ss_os_counts_and_ss_fraction_vs_pair_pt.png + .csv
// Usage : root -l -b -q 'plot_ss_fraction_vs_pair_pt.cxx+("tight")'
// =================================================================================================

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

namespace {

const std::string kDataDir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data";
const std::string kCtr     = "ctr50_80";   // requested Pb+Pb centrality bin (a nominal ctr_bin_edges bin)

// Muon working point (docs/muon_wp_registry.md): Tight = unsuffixed nominal file.
std::string WPFileSuffix(const std::string& wp)
{
    if (wp == "tight")  return "";
    if (wp == "medium") return "_medium_wp";
    throw std::runtime_error("plot_ss_fraction_vs_pair_pt: unknown WP '" + wp + "'");
}

// Pair-pT projection over ALL pair-eta bins of a (pT, eta) count histogram.
TH1D* ProjectPt(const std::string& path, const std::string& hname, const std::string& tag)
{
    TFile f(path.c_str(), "READ");
    if (f.IsZombie()) throw std::runtime_error("cannot open " + path);
    auto* h2 = dynamic_cast<TH2*>(f.Get(hname.c_str()));
    if (!h2) throw std::runtime_error("missing " + hname + " in " + path);
    // Data counts must be integer pair counts (Pb+Pb fills with the data "weight" column).
    if (std::fabs(h2->Integral(0, -1, 0, -1) - h2->GetEntries()) > 0.5)
        throw std::runtime_error(hname + " in " + path + ": integral != entries (non-unit weights)");
    TH1D* p = h2->ProjectionX(("p_" + tag).c_str(), 0, -1);
    p->SetDirectory(nullptr);
    return p;
}

void AssertSameAxis(const TH1* a, const TH1* b, const std::string& what)
{
    if (a->GetNbinsX() != b->GetNbinsX()) throw std::runtime_error("axis mismatch: " + what);
    for (int i = 1; i <= a->GetNbinsX() + 1; ++i)
        if (std::fabs(a->GetXaxis()->GetBinLowEdge(i) - b->GetXaxis()->GetBinLowEdge(i)) > 1e-9)
            throw std::runtime_error("axis edge mismatch: " + what);
}

TGraphAsymmErrors* Fraction(const TH1D* ss, const TH1D* os)
{
    auto* g = new TGraphAsymmErrors();
    for (int i = 1; i <= ss->GetNbinsX(); ++i) {
        const double n = ss->GetBinContent(i), N = n + os->GetBinContent(i);
        if (N <= 0) continue;
        const double f  = n / N;
        const double lo = TEfficiency::ClopperPearson(N, n, 0.683, false);
        const double hi = TEfficiency::ClopperPearson(N, n, 0.683, true);
        const int k = g->GetN();
        g->SetPoint(k, ss->GetBinCenter(i), 100. * f);
        g->SetPointError(k, ss->GetBinCenter(i) - ss->GetBinLowEdge(i),
                         ss->GetBinLowEdge(i + 1) - ss->GetBinCenter(i),
                         100. * (f - lo), 100. * (hi - f));
    }
    return g;
}

} // namespace

void plot_ss_fraction_vs_pair_pt(const std::string wp = "tight")
{
    const std::string wpsuf = WPFileSuffix(wp);

    // ---- pp24 ----
    const std::string pp_path = kDataDir + "/pp_2024/histograms_real_pairs_pp_2024_2mu4_nominal" + wpsuf + ".root";
    TH1D* pp_os = ProjectPt(pp_path, "h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts", "pp_os");
    TH1D* pp_ss = ProjectPt(pp_path, "h2d_ss_counts_pair_pt_pair_eta_binned_w_signal_cuts", "pp_ss");
    AssertSameAxis(pp_os, pp_ss, "pp OS vs SS");

    // ---- Pb+Pb 50-80 %, years summed as raw counts ----
    TH1D *aa_os = nullptr, *aa_ss = nullptr;
    for (int yr : {23, 24, 25, 26}) {
        const std::string ys = std::to_string(yr);
        const std::string path = kDataDir + "/pbpb_20" + ys + "/histograms_real_pairs_pbpb_20" + ys
                                 + "_single_mu4_no_trg_plots_nominal" + wpsuf + ".root";
        TH1D* os = ProjectPt(path, "h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + kCtr + "_counts", "aa_os" + ys);
        TH1D* ss = ProjectPt(path, "h2d_ss_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + kCtr + "_counts", "aa_ss" + ys);
        AssertSameAxis(pp_os, os, "pp vs PbPb" + ys + " OS");
        AssertSameAxis(pp_os, ss, "pp vs PbPb" + ys + " SS");
        std::cout << "PbPb" << ys << " " << kCtr << ": OS " << os->Integral(0, -1) << "  SS " << ss->Integral(0, -1) << "\n";
        if (!aa_os) { aa_os = os; aa_ss = ss; } else { aa_os->Add(os); aa_ss->Add(ss); }
    }

    // ---- CSV ----
    const std::string outdir = kDataDir + "/plots/single_b_analysis/ss_fraction";
    gSystem->mkdir(outdir.c_str(), true);
    const std::string base = outdir + "/ss_os_counts_and_ss_fraction_vs_pair_pt" + wpsuf;
    {
        std::ofstream csv(base + ".csv");
        csv << "pair_pt_lo,pair_pt_hi,pp_N_OS,pp_N_SS,pp_fSS_pct,PbPb5080_N_OS,PbPb5080_N_SS,PbPb5080_fSS_pct\n";
        csv << std::fixed;
        for (int i = 0; i <= pp_os->GetNbinsX() + 1; ++i) {   // incl. under/overflow rows
            auto pct = [](double s, double o) { return s + o > 0 ? 100. * s / (s + o) : NAN; };
            const double lo = i == 0 ? 0 : pp_os->GetBinLowEdge(i);
            const double hi = i > pp_os->GetNbinsX() ? INFINITY : pp_os->GetBinLowEdge(i + 1);
            csv << std::setprecision(2) << lo << "," << hi << "," << std::setprecision(0)
                << pp_os->GetBinContent(i) << "," << pp_ss->GetBinContent(i) << ","
                << std::setprecision(3) << pct(pp_ss->GetBinContent(i), pp_os->GetBinContent(i)) << ","
                << std::setprecision(0) << aa_os->GetBinContent(i) << "," << aa_ss->GetBinContent(i) << ","
                << std::setprecision(3) << pct(aa_ss->GetBinContent(i), aa_os->GetBinContent(i)) << "\n";
        }
    }

    // ---- plot ----
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TCanvas c("c", "", 1400, 1100);
    c.Divide(2, 2, 0.002, 0.002);

    const char* col_title[2] = {"#it{pp} 2024",
                                "Pb+Pb 2023#minus2026, 50#minus80 % centrality"};
    TH1D* os_h[2] = {pp_os, aa_os};
    TH1D* ss_h[2] = {pp_ss, aa_ss};
    const double xlo = pp_os->GetXaxis()->GetXmin(), xhi = pp_os->GetXaxis()->GetXmax();

    for (int col = 0; col < 2; ++col) {
        // top: raw counts. Log y: the counts fall over ~5 decades (power-law pT spectrum); log x
        // because the pair-pT axis is log-binned.
        c.cd(col + 1);
        gPad->SetLogx(); gPad->SetLogy();
        gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.08); gPad->SetBottomMargin(0.13);
        TH1D* os = os_h[col]; TH1D* ss = ss_h[col];
        os->SetLineColor(kBlue + 1); os->SetMarkerColor(kBlue + 1); os->SetMarkerStyle(20);
        ss->SetLineColor(kRed + 1);  ss->SetMarkerColor(kRed + 1);  ss->SetMarkerStyle(21);
        const double ymax = std::max(os->GetBinContent(os->GetMaximumBin()), ss->GetBinContent(ss->GetMaximumBin()));
        os->SetMinimum(0.5); os->SetMaximum(ymax * 8);
        os->GetXaxis()->SetTitle("p_{T}^{#mu#mu} [GeV]");
        os->GetYaxis()->SetTitle("Pairs (raw counts)");
        os->GetXaxis()->SetMoreLogLabels(); os->GetXaxis()->SetNoExponent();
        os->GetXaxis()->SetTitleSize(0.05); os->GetYaxis()->SetTitleSize(0.05);
        os->GetXaxis()->SetLabelSize(0.045); os->GetYaxis()->SetLabelSize(0.045);
        os->GetYaxis()->SetTitleOffset(1.3);
        os->Draw("E1");
        ss->Draw("E1 same");
        TLatex t; t.SetNDC(); t.SetTextSize(0.045);
        t.DrawLatex(0.16, 0.94, col_title[col]);
        auto* leg = new TLegend(0.62, 0.72, 0.94, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.045);
        leg->AddEntry(os, "opposite sign", "lep");
        leg->AddEntry(ss, "same sign", "lep");
        leg->Draw();

        // bottom: SS fraction, linear y.
        c.cd(col + 3);
        gPad->SetLogx();
        gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.04); gPad->SetBottomMargin(0.13);
        TGraphAsymmErrors* g = Fraction(ss, os);
        // Fixed range: the few-pair bins (N <= ~15) have Clopper-Pearson upper errors up to ~50 %,
        // which would squeeze every statistically meaningful point into the bottom of the pad; their
        // bars are clipped at the top (values in the CSV). Minimum slightly below 0 so the N_SS = 0
        // points stay visible above the frame line.
        auto* frame = new TH1D(Form("frame%d", col), "", 1, xlo, xhi);
        frame->SetMinimum(-0.5); frame->SetMaximum(20.);
        frame->GetXaxis()->SetTitle("p_{T}^{#mu#mu} [GeV]");
        frame->GetYaxis()->SetTitle("N_{SS} / (N_{SS} + N_{OS}) [%]");
        frame->GetXaxis()->SetMoreLogLabels(); frame->GetXaxis()->SetNoExponent();
        frame->GetXaxis()->SetTitleSize(0.05); frame->GetYaxis()->SetTitleSize(0.05);
        frame->GetXaxis()->SetLabelSize(0.045); frame->GetYaxis()->SetLabelSize(0.045);
        frame->GetYaxis()->SetTitleOffset(1.3);
        frame->Draw("AXIS");
        g->SetMarkerStyle(21); g->SetMarkerColor(kBlack); g->SetLineColor(kBlack);
        g->Draw("P same");
    }
    c.SaveAs((base + ".png").c_str());
    std::cout << "wrote " << base << ".png/.csv\n";
}
