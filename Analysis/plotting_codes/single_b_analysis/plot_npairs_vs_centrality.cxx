// N_pairs vs centrality (0-80%, 1% bins) for the PbPb running periods found on disk, combined.
// A period whose part files are not there yet (NTuple processing not run yet) is skipped with an
// [INFO] line and drops out of every title, legend entry and printed table -- TChain::Add only
// WARNS on a missing file, so without the probe the figure was silently overwritten with a title
// and a legend entry claiming a year that contributed nothing.
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

    struct YearInfo {
        int year; std::string dir; std::string tag; int nparts; int color; std::string label;
        std::vector<std::string> files;   // the part files that are actually on disk
    };
    // Per-year part counts = merged NTuple-processing part files on disk (2023 part5 / 2025 part7 =
    // Sep-2026 recovery skims; 2026 = 7).  Keep equal to file_batch_max in PbPbExtras.c;
    // pipelines/preflight_pbpb_year.sh <yr> cross-checks.  An UNDER-count silently reads a subset.
    std::vector<YearInfo> wanted = {
        {2023, base + "pbpb_2023/", "pbpb_2023", 5, kBlue+1,     "PbPb 2023", {}},
        {2024, base + "pbpb_2024/", "pbpb_2024", 2, kRed+1,      "PbPb 2024", {}},
        {2025, base + "pbpb_2025/", "pbpb_2025", 7, kGreen+2,    "PbPb 2025", {}},
        {2026, base + "pbpb_2026/", "pbpb_2026", 7, kMagenta+1,  "PbPb 2026", {}},
    };

    // Probe every expected part file. A year with no part on disk is dropped entirely; a year
    // with only some of its parts is kept, with the shortfall reported (a part whose NTuple
    // processing has not run yet is a normal transient state, not by itself an error).
    std::vector<YearInfo> years;
    for (auto& yr : wanted) {
        for (int p = 1; p <= yr.nparts; ++p) {
            std::string fname = yr.dir + "muon_pairs_" + yr.tag + "_part" + std::to_string(p) + suffix;
            if (!gSystem->AccessPathName(fname.c_str())) yr.files.push_back(fname);
        }
        if (yr.files.empty()) {
            printf("[INFO] %s: no input part file under %s -- skipping this year.\n",
                   yr.label.c_str(), yr.dir.c_str());
            continue;
        }
        printf("[INFO] %s: %zu of %d expected part files found.\n",
               yr.label.c_str(), yr.files.size(), yr.nparts);
        years.push_back(yr);
    }
    if (years.empty()) {
        printf("[FATAL] no PbPb input part file found for any year -- nothing to plot.\n");
        return;
    }

    // The histogram arrays hold one slot per SURVIVING year plus one extra "all years combined"
    // slot at index kComb. Both are derived from years.size() so that adding a year -- or
    // skipping one that is not on disk -- cannot leave a stale hard-coded 3/4 behind.
    const int kNY   = static_cast<int>(years.size());
    const int kComb = kNY;

    // "PbPb 2023+2024+2025", built from the years that actually contributed. Never a typed
    // year string: the figure must not be able to claim a year it did not read.
    std::string years_title = "PbPb ";
    for (int iy = 0; iy < kNY; ++iy) {
        if (iy) years_title += "+";
        years_title += std::to_string(years[iy].year);
    }

    const int nbins = 80;
    const double ctr_lo = 0, ctr_hi = 80;

    std::vector<TH1D*> h_os(kNY + 1), h_ss(kNY + 1);
    for (int i = 0; i <= kNY; ++i) {
        h_os[i] = new TH1D(Form("h_os_%d", i), "", nbins, ctr_lo, ctr_hi);
        h_ss[i] = new TH1D(Form("h_ss_%d", i), "", nbins, ctr_lo, ctr_hi);
    }

    for (int iy = 0; iy < kNY; ++iy) {
        auto& yr = years[iy];
        TChain chain_os("muon_pair_tree_sign2");  // OS
        TChain chain_ss("muon_pair_tree_sign1");  // SS
        for (const auto& fname : yr.files) {
            chain_os.Add(fname.c_str());
            chain_ss.Add(fname.c_str());
        }
        printf("%s: OS entries = %lld, SS entries = %lld\n",
               yr.label.c_str(), chain_os.GetEntries(), chain_ss.GetEntries());

        fill_chain(chain_os, h_os[iy]);
        fill_chain(chain_ss, h_ss[iy]);

        h_os[kComb]->Add(h_os[iy]);
        h_ss[kComb]->Add(h_ss[iy]);
    }

    // Print pair counts
    printf("\n--- Pair counts in 0-80%% centrality ---\n");
    for (int iy = 0; iy < kNY; ++iy)
        printf("%-15s  OS = %.0f   SS = %.0f\n",
               years[iy].label.c_str(), h_os[iy]->Integral(), h_ss[iy]->Integral());
    printf("%-15s  OS = %.0f   SS = %.0f\n", "Combined",
           h_os[kComb]->Integral(), h_ss[kComb]->Integral());

    // --- Canvas: 2 rows × 2 cols ---
    TCanvas* c = new TCanvas("c", "", 1200, 900);
    c->Divide(2, 2);

    auto draw_log = [&](int pad, std::vector<TH1D*>& harr, const char* label) {
        c->cd(pad);
        gPad->SetLogy();
        TH1D* hc = harr[kComb];
        hc->SetLineColor(kBlack); hc->SetLineWidth(2);
        hc->GetXaxis()->SetTitle("Centrality (%)");
        hc->GetYaxis()->SetTitle(Form("N_{pairs} (%s) / 1%%", label));
        hc->SetTitle(Form("%s, %s", years_title.c_str(), label));
        hc->SetMinimum(0.5); hc->SetMaximum(hc->GetMaximum() * 5);
        hc->Draw("hist");
        for (int iy = 0; iy < kNY; ++iy) {
            harr[iy]->SetLineColor(years[iy].color);
            harr[iy]->SetLineWidth(1);
            harr[iy]->SetMinimum(0.5);
            harr[iy]->Draw("hist same");
        }
        // y-range grows with the number of year entries (combined + one per year)
        TLegend* leg = new TLegend(0.55, 0.88 - 0.065 * (kNY + 1), 0.88, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0);
        leg->AddEntry(hc, "Combined", "l");
        for (int iy = 0; iy < kNY; ++iy) leg->AddEntry(harr[iy], years[iy].label.c_str(), "l");
        leg->Draw();
    };

    draw_log(1, h_os, "OS");
    draw_log(2, h_ss, "SS");

    // Pad 3: combined OS, linear (clone to avoid log-scale contamination from pad 1)
    c->cd(3);
    gPad->SetLeftMargin(0.15);   // the y-axis title was clipped by the default margin
    TH1D* h_os_lin = (TH1D*)h_os[kComb]->Clone("h_os_lin");
    h_os_lin->SetTitle((years_title + " combined, OS (linear)").c_str());
    // NB: TH1::GetMaximum() returns the STORED fMaximum when one has been set, not the
    // largest bin content -- and this clone inherited the 5x log-scale maximum that
    // draw_log() set on h_os[kComb].  Using it here just rescaled that stored value, so
    // the linear pad was drawn on a ~6e4 axis for a ~1e4 peak.  Read the bin content.
    h_os_lin->SetMinimum(0);
    h_os_lin->SetMaximum(h_os_lin->GetBinContent(h_os_lin->GetMaximumBin()) * 1.2);
    h_os_lin->Draw("hist");

    // Pad 4: OS/SS
    c->cd(4);
    gPad->SetLeftMargin(0.15);   // the y-axis title was clipped by the default margin
    TH1D* h_ratio = (TH1D*)h_os[kComb]->Clone("h_ratio");
    h_ratio->Divide(h_ss[kComb]);
    h_ratio->SetTitle("OS / SS (combined)");
    h_ratio->GetYaxis()->SetTitle("OS / SS");
    // The clone inherits the axis range draw_log() set on h_os[kComb] (max = 5x the OS
    // peak, i.e. ~5e4), so the ratio -- an O(1) quantity -- was drawn far below the
    // visible area and pad 4 rendered EMPTY.  (Pre-existing since the panel was written;
    // unrelated to the year handling.)  GetMaximum() cannot be used to fix it: it returns
    // the STORED fMaximum once SetMaximum has been called, so it would just hand back
    // that same 5e4.  Take the largest bin content instead.
    h_ratio->SetMinimum(0.);
    h_ratio->SetMaximum(h_ratio->GetBinContent(h_ratio->GetMaximumBin()) * 1.2);
    h_ratio->Draw("hist");

    const std::string outdir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/";
    gSystem->MakeDirectory(outdir.c_str());
    c->SaveAs((outdir + "npairs_vs_centrality_pbpb_combined.png").c_str());
    printf("Saved to %s\n", outdir.c_str());
}
