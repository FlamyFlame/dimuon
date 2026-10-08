// Abnormal-run check: (# signal-region events) / (integrated luminosity) per run, one PNG per dataset.
//   Numerator  : OS pairs of the NOMINAL muon-pair NTuple passing the crossx signal selection, Tight WP
//                (PbPb: NTuple stage already applies the 5-cut event selection; + centrality acceptance
//                [0,80) as in FillHistogramsCrossx). RAW counts -- no efficiency weights.
//                Counted as distinct events (run, LB, BCID, entry); differs from the pair count by <0.3 %.
//   Denominator: 'Prescale Corrected' column of the lumi table of the analysis trigger
//                (pp24: HLT_2mu4_L12MU3V in pb^-1 ; PbPb: HLT_mu4, file unit ub^-1 -> nb^-1).
// Runs excluded at event level (PbPbBadRuns) and runs with zero luminosity are not drawn.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "NTupleProcessingCode/PbPbEventSelConfig.h"
#include "plotting_codes/data_qa/per_run_signal_yield_table.cxx"

namespace {
struct DS { std::string tag, label, file, lumi_csv, lumi_unit; bool pbpb; int year; double lumi_scale; };
const std::string DATA = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
const std::string LUMI = "/gpfs/mnt/atlasgpfs01/usatlas/workarea/yuhanguo/dimuon_codes/IntNotes/data/luminosity/";
}

static std::map<int,double> ReadLumi(const std::string& path) {
    std::map<int,double> m; std::ifstream f(path); std::string line; std::getline(f, line);
    while (std::getline(f, line)) {
        std::stringstream ss(line); std::string c; std::vector<std::string> v;
        while (std::getline(ss, c, ',')) v.push_back(c);
        if (v.size() < 7 || v[0].find_first_not_of("0123456789 ") != std::string::npos) continue; // skip Total row
        m[std::stoi(v[0])] = std::stod(v[6]);                                                    // Prescale Corrected
    }
    return m;
}

static void PlotOne(const DS& d, const std::string& outdir) {
    per_run_signal_yield_table(d.tag, d.file, d.pbpb, outdir);
    const auto lumi = ReadLumi(LUMI + d.lumi_csv);
    const auto bad  = d.pbpb ? PbPbBadRuns(d.year) : std::unordered_set<int>{};
    std::map<int,double> N;
    { std::ifstream f(outdir + "/" + d.tag + "_per_run.csv"); std::string l; std::getline(f, l);
      while (std::getline(f, l)) { int r; long np, ne; double sw; char c;
          std::stringstream ss(l); ss >> r >> c >> np >> c >> ne >> c >> sw; N[r] = ne; } }
    std::vector<int> runs; double sN = 0, sL = 0;
    for (auto& kv : lumi) { if (kv.second <= 0 || bad.count(kv.first)) continue; runs.push_back(kv.first);
                            sN += N.count(kv.first) ? N[kv.first] : 0; sL += kv.second * d.lumi_scale; }
    for (auto& kv : N) if (!lumi.count(kv.first)) throw std::runtime_error("run " + std::to_string(kv.first) + " in data but not in lumi table");
    std::sort(runs.begin(), runs.end());
    const double mean = sN / sL; const int nr = runs.size();

    TH1D* fr = new TH1D(("fr_" + d.tag).c_str(), "", nr, 0, nr);
    std::vector<double> y(nr), ey(nr);
    std::ofstream out(outdir + "/" + d.tag + "_yield_per_lumi.csv");
    out << "run,N_events,lumi_" << d.lumi_unit << ",N_over_lumi,err,pull_vs_mean\n";
    double ymax = 0;
    for (int i = 0; i < nr; ++i) {
        const int r = runs[i]; const double L = lumi.at(r) * d.lumi_scale, n = N.count(r) ? N[r] : 0;
        y[i] = n / L; ey[i] = std::sqrt(std::max(n, 1.0)) / L;
        out << r << "," << n << "," << L << "," << y[i] << "," << ey[i] << "," << (n - mean * L) / std::sqrt(mean * L) << "\n";
        fr->GetXaxis()->SetBinLabel(i + 1, std::to_string(r).c_str());
        ymax = std::max(ymax, y[i] + ey[i]);
    }
    gStyle->SetOptStat(0); TGaxis::SetMaxDigits(6);
    TCanvas* c = new TCanvas(("c_" + d.tag).c_str(), "", 1400, 750);
    c->SetBottomMargin(0.16); c->SetLeftMargin(0.09); c->SetRightMargin(0.03); c->SetTopMargin(0.07);
    fr->SetMinimum(0); fr->SetMaximum(1.35 * ymax);
    fr->GetXaxis()->SetLabelSize(nr > 40 ? 0.022 : 0.03); fr->GetXaxis()->LabelsOption("v"); fr->GetXaxis()->SetTickLength(0);
    fr->GetXaxis()->SetTitle("Run number"); fr->GetXaxis()->SetTitleOffset(2.2); fr->GetXaxis()->SetTitleSize(0.04);
    fr->GetYaxis()->SetTitle(Form("N_{events}^{signal} / L  [%s^{-1}]", d.lumi_unit.c_str()));
    fr->GetYaxis()->SetTitleSize(0.04); fr->GetYaxis()->SetTitleOffset(1.1);
    fr->Draw("AXIS");
    TGraphAsymmErrors* g = new TGraphAsymmErrors(nr);
    for (int i = 0; i < nr; ++i) { g->SetPoint(i, i + 0.5, y[i]); g->SetPointError(i, 0, 0, ey[i], ey[i]); }
    g->SetMarkerStyle(20); g->SetMarkerSize(1.0); g->SetLineWidth(2);
    TLine* ln = new TLine(0, mean, nr, mean); ln->SetLineStyle(2); ln->SetLineColor(kRed + 1); ln->SetLineWidth(2);
    ln->Draw(); g->Draw("P SAME");
    TLegend* lg = new TLegend(0.60, 0.78, 0.96, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0);
    lg->AddEntry(g, "per run (stat. error)", "pe");
    lg->AddEntry(ln, Form("dataset average %.4g %s^{-1}", mean, d.lumi_unit.c_str()), "l"); lg->Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.04); t.SetTextFont(72); t.DrawLatex(0.10, 0.94, "ATLAS");
    t.SetTextFont(42); t.DrawLatex(0.17, 0.94, ("Internal   " + d.label).c_str());
    c->SaveAs((outdir + "/" + d.tag + "_signal_yield_per_lumi_per_run.png").c_str());
}

void plot_signal_yield_per_lumi_per_run() {
    const std::string outdir = DATA + "plots/data_qa/signal_yield_per_lumi_per_run";
    gSystem->mkdir(outdir.c_str(), true);
    std::vector<DS> ds = {
        {"pp24",  "pp 2024, HLT_2mu4_L12MU3V, #sqrt{s} = 5.36 TeV",  DATA + "pp_2024/muon_pairs_pp_2024_2mu4_mindR_0_02.root",              "pp_2024/lumitable_pp_24_HLT_2mu4_L12MU3V.csv", "pb",  false, 2024, 1.0},
        {"pbpb23","Pb+Pb 2023, HLT_mu4, #sqrt{s_{NN}} = 5.36 TeV",   DATA + "pbpb_2023/muon_pairs_pbpb_2023_single_mu4_mindR_0_02.root",    "pbpb_2023/lumitable_pbpb_23_HLT_mu4.csv",      "nb",  true, 2023, 1e-3},
        {"pbpb24","Pb+Pb 2024, HLT_mu4, #sqrt{s_{NN}} = 5.36 TeV",   DATA + "pbpb_2024/muon_pairs_pbpb_2024_single_mu4_mindR_0_02.root",    "pbpb_2024/lumitable_pbpb_24_HLT_mu4.csv",      "nb",  true, 2024, 1e-3},
        {"pbpb25","Pb+Pb 2025, HLT_mu4, #sqrt{s_{NN}} = 5.36 TeV",   DATA + "pbpb_2025/muon_pairs_pbpb_2025_single_mu4_mindR_0_02.root",    "pbpb_2025/lumitable_pbpb_25_HLT_mu4.csv",      "nb",  true, 2025, 1e-3},
        {"pbpb26","Pb+Pb 2026, HLT_mu4, #sqrt{s_{NN}} = 5.36 TeV",   DATA + "pbpb_2026/muon_pairs_pbpb_2026_single_mu4_mindR_0_02.root",    "pbpb_2026/lumitable_pbpb_26_HLT_mu4.csv",      "nb",  true, 2026, 1e-3},
    };
    for (auto& d : ds) PlotOne(d, outdir);
}
