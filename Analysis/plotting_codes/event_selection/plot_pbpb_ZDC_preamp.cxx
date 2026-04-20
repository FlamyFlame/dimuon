// plot_pbpb_ZDC_preamp.cxx
// Overlay ZDC per-arm presample amplitude (sum of 4 modules) for A and C sides.
// Reads directly from the pair ntuple (sign1 + sign2 trees).
//
// Usage:
//   .L plot_pbpb_ZDC_preamp.cxx+
//   plot_pbpb_ZDC_preamp(24)

#include <string>
#include <stdexcept>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TLegend.h"

void plot_pbpb_ZDC_preamp(int run_year = 24) {
    run_year %= 2000;
    const std::string yr = std::to_string(run_year);

    const std::string infile =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20" + yr +
        "/muon_pairs_pbpb_20" + yr + "_single_mu4_mindR_0_02_res_cut_v2.root";
    const std::string out_dir =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/event_selection";

    if (gSystem->AccessPathName(infile.c_str()))
        throw std::runtime_error("Input file not found: " + infile);
    TFile* f = TFile::Open(infile.c_str(), "READ");
    if (!f || f->IsZombie())
        throw std::runtime_error("Cannot open: " + infile);

    gSystem->mkdir(out_dir.c_str(), true);

    // 150 bins, -1000 to 6500
    TH1D* hA = new TH1D("hA", ";ZDC preamp sum [ADC counts];Pairs", 150, -1000, 6500);
    TH1D* hC = new TH1D("hC", ";ZDC preamp sum [ADC counts];Pairs", 150, -1000, 6500);
    hA->SetDirectory(nullptr);
    hC->SetDirectory(nullptr);

    for (const char* tname : {"muon_pair_tree_sign1", "muon_pair_tree_sign2"}) {
        TTree* t = (TTree*)f->Get(tname);
        if (!t) { std::cerr << tname << " not found, skipping" << std::endl; continue; }
        t->SetMakeClass(1);
        Float_t pa = -999.f, pc = -999.f;
        t->SetBranchAddress("MuonPairObj.ZDC_preamp_A", &pa);
        t->SetBranchAddress("MuonPairObj.ZDC_preamp_C", &pc);
        const Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            hA->Fill(pa);
            hC->Fill(pc);
        }
    }
    f->Close();

    std::cout << "preamp_A entries=" << hA->GetEntries()
              << "  preamp_C entries=" << hC->GetEntries() << std::endl;

    // No normalisation — raw pair counts

    hA->SetLineColor(kBlue);
    hA->SetLineWidth(2);
    hA->SetMarkerColor(kBlue);
    hA->SetMarkerStyle(20);
    hA->SetMarkerSize(0.5);

    hC->SetLineColor(kRed);
    hC->SetLineWidth(2);
    hC->SetMarkerColor(kRed);
    hC->SetMarkerStyle(20);
    hC->SetMarkerSize(0.5);

    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas("c_preamp", "", 800, 620);
    c->SetLeftMargin(0.13);
    c->SetBottomMargin(0.13);
    c->SetLogy();

    double ymax = 3.0 * std::max(hA->GetMaximum(), hC->GetMaximum());
    hA->GetYaxis()->SetRangeUser(0.5, ymax);
    hA->Draw("E");
    hC->Draw("E SAME");

    TLegend* leg = new TLegend(0.55, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.036);
    leg->AddEntry(hA, "ZDC preamp A-side", "lp");
    leg->AddEntry(hC, "ZDC preamp C-side", "lp");
    leg->Draw();

    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.036);
    tl.DrawLatex(0.15, 0.87, Form("Pb+Pb 20%s data", yr.c_str()));

    const std::string outpath = out_dir + "/ZDC_preamp_A_vs_C_pbpb_20" + yr + ".png";
    c->SaveAs(outpath.c_str());
    std::cout << "Saved: " << outpath << std::endl;
    delete c;
    delete hA;
    delete hC;
}
