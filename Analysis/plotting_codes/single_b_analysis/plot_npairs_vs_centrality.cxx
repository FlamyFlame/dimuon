// N_pairs vs centrality (0-80%, 1% bins) for PbPb 23+24+25 combined
// Uses SetMakeClass(1) to read avg_centrality without needing a compiled dictionary
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include <vector>
#include <string>
#include <iostream>

static void fill_chain(TChain& chain, TH1D* h) {
    chain.SetMakeClass(1);
    Int_t avg_centrality = -1;
    // with SetMakeClass(1), split branches are flattened: leaf name = "avg_centrality"
    int rc = chain.SetBranchAddress("avg_centrality", &avg_centrality);
    if (rc < 0) {
        std::cerr << "ERROR: SetBranchAddress returned " << rc << " for avg_centrality\n";
        return;
    }
    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("avg_centrality", 1);

    Long64_t n = chain.GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        chain.GetEntry(i);
        if (avg_centrality >= 0 && avg_centrality < 80)
            h->Fill(avg_centrality + 0.5);  // centre in the 1% bin
    }
}

void plot_npairs_vs_centrality() {
    gStyle->SetOptStat(0);
    gROOT->SetBatch(true);

    const std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
    const std::string suffix = "_single_mu4_mindR_0_01_res_cut_v2.root";

    struct YearInfo { std::string dir; std::string tag; int nparts; int color; std::string label; };
    std::vector<YearInfo> years = {
        {base + "pbpb_2023/", "pbpb_2023", 4, kBlue+1,  "PbPb 2023"},
        {base + "pbpb_2024/", "pbpb_2024", 2, kRed+1,   "PbPb 2024"},
        {base + "pbpb_2025/", "pbpb_2025", 6, kGreen+2, "PbPb 2025"},
    };

    const int nbins = 80;
    const double ctr_lo = 0, ctr_hi = 80;

    TH1D* h_os[4];
    TH1D* h_ss[4];
    for (int i = 0; i < 4; ++i) {
        h_os[i] = new TH1D(Form("h_os_%d", i), "", nbins, ctr_lo, ctr_hi);
        h_ss[i] = new TH1D(Form("h_ss_%d", i), "", nbins, ctr_lo, ctr_hi);
    }

    for (int iy = 0; iy < (int)years.size(); ++iy) {
        auto& yr = years[iy];
        TChain chain_os("muon_pair_tree_sign2");  // OS
        TChain chain_ss("muon_pair_tree_sign1");  // SS
        for (int p = 1; p <= yr.nparts; ++p) {
            std::string fname = yr.dir + "muon_pairs_" + yr.tag + "_part" + std::to_string(p) + suffix;
            chain_os.Add(fname.c_str());
            chain_ss.Add(fname.c_str());
        }
        printf("%s: OS entries = %lld, SS entries = %lld\n",
               yr.label.c_str(), chain_os.GetEntries(), chain_ss.GetEntries());

        fill_chain(chain_os, h_os[iy]);
        fill_chain(chain_ss, h_ss[iy]);

        h_os[3]->Add(h_os[iy]);
        h_ss[3]->Add(h_ss[iy]);
    }

    // Print pair counts
    printf("\n--- Pair counts in 0-80%% centrality ---\n");
    for (int iy = 0; iy < 3; ++iy)
        printf("%-15s  OS = %.0f   SS = %.0f\n",
               years[iy].label.c_str(), h_os[iy]->Integral(), h_ss[iy]->Integral());
    printf("%-15s  OS = %.0f   SS = %.0f\n", "Combined", h_os[3]->Integral(), h_ss[3]->Integral());

    // --- Canvas: 2 rows × 2 cols ---
    TCanvas* c = new TCanvas("c", "", 1200, 900);
    c->Divide(2, 2);

    auto draw_log = [&](int pad, TH1D** harr, const char* label) {
        c->cd(pad);
        gPad->SetLogy();
        TH1D* hc = harr[3];
        hc->SetLineColor(kBlack); hc->SetLineWidth(2);
        hc->GetXaxis()->SetTitle("Centrality (%)");
        hc->GetYaxis()->SetTitle(Form("N_{pairs} (%s) / 1%%", label));
        hc->SetTitle(Form("PbPb 23+24+25, %s", label));
        hc->SetMinimum(0.5); hc->SetMaximum(hc->GetMaximum() * 5);
        hc->Draw("hist");
        for (int iy = 0; iy < 3; ++iy) {
            harr[iy]->SetLineColor(years[iy].color);
            harr[iy]->SetLineWidth(1);
            harr[iy]->SetMinimum(0.5);
            harr[iy]->Draw("hist same");
        }
        TLegend* leg = new TLegend(0.55, 0.62, 0.88, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->AddEntry(hc, "Combined", "l");
        for (int iy = 0; iy < 3; ++iy) leg->AddEntry(harr[iy], years[iy].label.c_str(), "l");
        leg->Draw();
    };

    draw_log(1, h_os, "OS");
    draw_log(2, h_ss, "SS");

    // Pad 3: combined OS, linear (clone to avoid log-scale contamination from pad 1)
    c->cd(3);
    TH1D* h_os_lin = (TH1D*)h_os[3]->Clone("h_os_lin");
    h_os_lin->SetTitle("PbPb 23+24+25 combined, OS (linear)");
    h_os_lin->SetMinimum(0);
    h_os_lin->SetMaximum(h_os_lin->GetMaximum() * 1.2);
    h_os_lin->Draw("hist");

    // Pad 4: OS/SS
    c->cd(4);
    TH1D* h_ratio = (TH1D*)h_os[3]->Clone("h_ratio");
    h_ratio->Divide(h_ss[3]);
    h_ratio->SetTitle("OS / SS (combined)");
    h_ratio->GetYaxis()->SetTitle("OS / SS");
    h_ratio->SetMinimum(0);
    h_ratio->Draw("hist");

    const std::string outdir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/";
    gSystem->MakeDirectory(outdir.c_str());
    c->SaveAs((outdir + "npairs_vs_centrality_pbpb_combined.png").c_str());
    printf("Saved to %s\n", outdir.c_str());
}
