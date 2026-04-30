// plot_zdc_fcal_per_run_pbpb23.cxx
//
// ZDC total energy vs FCal ET 2D, one pad per run, all 60 pbpb2023 runs.
// Axes: x = FCal ET^{A+C} in [-0.5, 5.5] TeV, y = ZDC E_total in [0, 400] TeV.
// Color scale: log z.  Selection: b_HLT_mu4_L1MU3V only.
// Output: 3 PNG files, 20 runs each (5 col × 4 row per file).
//
// Usage (run from Analysis/):
//   source /usatlas/u/yuhanguo/setup.sh
//   root -l -b -q 'plotting_codes/event_selection/plot_zdc_fcal_per_run_pbpb23.cxx+'

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include "TChain.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TStyle.h"

static const std::string kBase    = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
static const std::string kPlotDir = kBase + "plots/single_b_analysis/event_selection/pbpb_2023/run_dependence/";

static const int kNCols = 5, kNRows = 4, kPerFile = kNCols * kNRows;  // 20 runs per PNG

static void DrawRunGrid(const std::vector<int>& runs, int start,
                        const std::map<int, TH2D*>& hmap,
                        const std::map<int, long long>& nmap,
                        const std::string& outPath) {
    int n = std::min(kPerFile, (int)runs.size() - start);
    TCanvas* cv = new TCanvas(Form("cv_gr%d", start), "", 1600, 1280);
    cv->SetFillColor(0);

    for (int idx = 0; idx < n; ++idx) {
        int run = runs[start + idx];
        int row = idx / kNCols;
        int col = idx % kNCols;

        double xLo = col              * (1.0 / kNCols);
        double xHi = (col + 1)        * (1.0 / kNCols);
        double yLo = (kNRows-row-1)   * (1.0 / kNRows);
        double yHi = (kNRows-row)     * (1.0 / kNRows);

        TPad* pad = new TPad(Form("pad%d_%d", start, idx), "", xLo, yLo, xHi, yHi);
        pad->SetLeftMargin(0.13);
        pad->SetRightMargin(0.17);
        pad->SetTopMargin(0.08);
        pad->SetBottomMargin(0.18);
        pad->SetLogz();
        pad->Draw();
        pad->cd();

        TH2D* h = hmap.at(run);
        h->SetContour(99);
        h->GetXaxis()->SetTitleSize(0.065);
        h->GetXaxis()->SetLabelSize(0.050);
        h->GetYaxis()->SetTitleSize(0.065);
        h->GetYaxis()->SetLabelSize(0.050);
        h->GetYaxis()->SetTitleOffset(0.90);
        h->GetZaxis()->SetLabelSize(0.048);
        h->SetMinimum(0.5);
        h->GetXaxis()->SetRangeUser(-0.5, 5.5);
        h->GetYaxis()->SetRangeUser(0., 400.);
        h->Draw("COLZ");

        long long nev = nmap.count(run) ? nmap.at(run) : 0LL;
        TLatex tl; tl.SetNDC(); tl.SetTextAlign(11);
        tl.SetTextSize(0.070);
        tl.DrawLatex(0.15, 0.90, Form("Run %d", run));
        tl.SetTextSize(0.060);
        tl.DrawLatex(0.15, 0.81, Form("N = %lld", nev));

        cv->cd();
    }
    cv->SaveAs(outPath.c_str());
    std::cout << "Saved: " << outPath << std::endl;
    delete cv;
}

void plot_zdc_fcal_per_run_pbpb23() {
    gStyle->SetOptStat(0);
    gSystem->mkdir(kPlotDir.c_str(), true);

    // ---- Build TChain ----
    TChain chain("HeavyIonD3PD", "HeavyIonD3PD");
    for (int p = 1; p <= 4; ++p) {
        std::string f = kBase + "pbpb_2023/data_pbpb23_part" + std::to_string(p) + ".root";
        if (!gSystem->AccessPathName(f.c_str())) chain.Add(f.c_str());
        else std::cerr << "Skipping: " << f << std::endl;
    }
    chain.SetMakeClass(1);
    chain.SetBranchStatus("*", 0);

    Int_t   b_RunNum = 0;
    Bool_t  b_HLT   = false;
    Float_t FCal_P  = 0.f, FCal_N = 0.f;
    Float_t zdcE[2]{};

    chain.SetBranchStatus("RunNumber",         1);
    chain.SetBranchStatus("b_HLT_mu4_L1MU3V",  1);
    chain.SetBranchStatus("FCal_Et_P",         1);
    chain.SetBranchStatus("FCal_Et_N",         1);
    chain.SetBranchStatus("zdc_ZdcEnergy",     1);
    chain.SetBranchAddress("RunNumber",         &b_RunNum);
    chain.SetBranchAddress("b_HLT_mu4_L1MU3V", &b_HLT);
    chain.SetBranchAddress("FCal_Et_P",         &FCal_P);
    chain.SetBranchAddress("FCal_Et_N",         &FCal_N);
    chain.SetBranchAddress("zdc_ZdcEnergy",     zdcE);

    std::map<int, TH2D*>     hmap;
    std::map<int, long long> nmap;

    const Long64_t nTot = chain.GetEntries();
    std::cout << "Scanning " << nTot << " events..." << std::flush;
    for (Long64_t i = 0; i < nTot; ++i) {
        chain.GetEntry(i);
        if (!b_HLT) continue;

        int run = b_RunNum;
        if (!hmap.count(run)) {
            hmap[run] = new TH2D(Form("hEF_%d", run),
                                 ";FCal E_{T}^{A+C} [TeV];ZDC E_{total} [TeV]",
                                 120, -0.5, 5.5, 100, 0., 400.);
            hmap[run]->SetDirectory(nullptr);
            nmap[run] = 0LL;
        }
        const float fcal   = (FCal_P + FCal_N) * 1e-6f;
        const float zdcTot = (zdcE[0] + zdcE[1]) / 1000.f;
        hmap[run]->Fill(fcal, zdcTot);
        ++nmap[run];
    }
    std::cout << "  done. " << hmap.size() << " runs found.\n";

    // Sorted run list
    std::vector<int> runs;
    runs.reserve(hmap.size());
    for (auto& p : hmap) runs.push_back(p.first);
    std::sort(runs.begin(), runs.end());

    for (int file = 0; file < 3; ++file) {
        int start = file * kPerFile;
        if (start >= (int)runs.size()) break;
        std::string outPath = kPlotDir + Form("zdc_fcal_per_run_pbpb23_part%d.png", file + 1);
        DrawRunGrid(runs, start, hmap, nmap, outPath);
    }

    for (auto& p : hmap) delete p.second;
}
