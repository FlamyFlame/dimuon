// plot_dr_vs_pair_pt_diagnostic.cxx
//
// Diagnostic: 2D map of dimuon DeltaR (zoomed 0-1, so the small-DeltaR region
// BELOW 0.05 is visible) vs pair pT (log x-axis), applying ALL single-b
// signal-selection cuts EXCEPT the DeltaR cut. Run 3 single-b dimuon analysis.
//
// Self-contained standalone macro (no pipeline class touched). Reads the
// muon_pairs trees directly, structured after
//   plotting_codes/single_b_analysis/plot_sig_accept_cutflow_above_60GeV.cxx
//
// Three plots (each TH2 colz, x = pair pT log axis w/ analysis pT_bins_120
// binning, y = DeltaR uniform 0-1, 50 bins; PNG only):
//   1. pythia_truth_dr_vs_pair_pt.png  (Pythia evgen single-b, opposite-sign, weighted)
//   2. pp24_data_dr_vs_pair_pt.png     (pp24 data, raw counts)
//   3. pbpb_data_dr_vs_pair_pt.png     (PbPb 23+24+25 combined, raw counts)
//
// A horizontal dashed red line is drawn at DeltaR = 0.05 (the removed cut).
//
// Run (ACLiC; needed for struct-member access in JIT filters):
//   root -l -b -q 'plotting_codes/single_b_analysis/plot_dr_vs_pair_pt_diagnostic.cxx+'

#include "../../MuonObjectsParamsAndHelpers/MuonPairPythia.h"  // truth pair struct (pulls in MuonPairPbPb/Reco)
#include "../../MuonObjectsParamsAndHelpers/MuonPairPbPb.h"    // PbPb data pair struct
#include "../../MuonObjectsParamsAndHelpers/MuonPairReco.h"    // PP data pair struct (MuonPairPP)

#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH2D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// Replicate ParamsSet::fillLogBinningArray + pT_bins_120 (15 log bins, 8 -> 120 GeV).
static std::vector<double> MakeLogBins(int nBins, double low, double high) {
    std::vector<double> bins;
    const double logLow  = std::log10(low);
    const double logHigh = std::log10(high);
    const double logStep = (logHigh - logLow) / nBins;
    for (int i = 0; i <= nBins; ++i) bins.push_back(std::pow(10.0, logLow + i * logStep));
    return bins;
}

static const double kDrCut = 0.05;  // the removed DeltaR cut value

// Draw one TH2 colz with log x, dashed red line at DeltaR = 0.05, save PNG.
static void DrawAndSave(TH2D& h, const std::string& out_path, const std::string& header)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TCanvas c(("c_" + std::string(h.GetName())).c_str(), "", 950, 720);
    c.SetLogx();
    c.SetRightMargin(0.15);
    c.SetLeftMargin(0.12);

    // Large dynamic range in raw counts / weights -> log z.
    if (h.GetMaximum() > 0) c.SetLogz();

    h.GetXaxis()->SetTitle("p_{T}^{#mu#mu} [GeV]");
    h.GetYaxis()->SetTitle("#DeltaR(#mu_{1},#mu_{2})");
    h.GetXaxis()->SetMoreLogLabels();
    h.GetXaxis()->SetNoExponent();
    h.Draw("colz");

    // Horizontal dashed red line at the removed DeltaR cut.
    const double xlo = h.GetXaxis()->GetXmin();
    const double xhi = h.GetXaxis()->GetXmax();
    TLine line(xlo, kDrCut, xhi, kDrCut);
    line.SetLineColor(kRed);
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.Draw();

    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.032);
    lat.DrawLatex(0.12, 0.945, header.c_str());
    lat.DrawLatex(0.12, 0.905, "all single-b cuts except #DeltaR;  red dashed = #DeltaR=0.05");

    c.SaveAs(out_path.c_str());
    std::cout << "[INFO] Saved: " << out_path << "\n";
}

void plot_dr_vs_pair_pt_diagnostic()
{
    ROOT::EnableImplicitMT();

    // ---- axis binning ----
    const std::vector<double> xbins = MakeLogBins(15, 8.0, 120.0);  // pT_bins_120
    const int    nx = static_cast<int>(xbins.size()) - 1;
    const int    ny = 50;
    const double ylo = 0.0, yhi = 1.0;

    const std::string out_dir =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis/dr_vs_pair_pt_diagnostic";
    gSystem->mkdir(out_dir.c_str(), true);

    // ---- inputs ----
    const std::string truth_file =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/"
        "muon_pairs_pythia_5p36TeV_no_data_resonance_cuts.root";

    // pp24 nominal crossx input: trigger_mode=3 (2mu4), mindR_0_02, combined/hadded.
    const std::string pp24_file =
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/"
        "muon_pairs_pp_2024_2mu4_mindR_0_02.root";

    // PbPb nominal crossx inputs: trigger_mode=1 (single_mu4), mindR_0_02, res_cut_v2, per-year hadded.
    const std::vector<std::string> pbpb_files = {
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/"
        "muon_pairs_pbpb_2023_single_mu4_mindR_0_02_res_cut_v2.root",
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2024/"
        "muon_pairs_pbpb_2024_single_mu4_mindR_0_02_res_cut_v2.root",
        "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2025/"
        "muon_pairs_pbpb_2025_single_mu4_mindR_0_02_res_cut_v2.root",
    };

    auto check = [](const std::string& f) {
        if (gSystem->AccessPathName(f.c_str())) {
            std::cerr << "[FATAL] missing input file: " << f << "\n";
            gSystem->Exit(1);
        }
    };
    check(truth_file);
    check(pp24_file);
    for (auto& f : pbpb_files) check(f);

    // Signal cuts EXCEPT DeltaR.
    const std::string truth_cuts =
        "from_same_b && truth_minv > 1.08 && truth_minv < 2.9 && truth_pair_pt > 8 "
        "&& m1.truth_charge * m1.truth_eta < 2.2 && m2.truth_charge * m2.truth_eta < 2.2";
    const std::string data_cuts =
        "minv > 1.08 && minv < 2.9 && pair_pt > 8 "
        "&& m1.charge * m1.eta < 2.2 && m2.charge * m2.eta < 2.2";

    // ======================= 1. TRUTH (Pythia single-b) =======================
    {
        std::cout << "\n[INFO] Pythia truth single-b...\n";
        auto df = ROOT::RDataFrame("muon_pair_tree_sign2", std::vector<std::string>{truth_file})
                      .Filter(truth_cuts);
        ROOT::RDF::TH2DModel model("h_truth", "Pythia truth single-b;p_{T}^{#mu#mu} [GeV];#DeltaR",
                                   nx, xbins.data(), ny, ylo, yhi);
        auto h = df.Histo2D(model, "truth_pair_pt", "truth_dr", "weight");
        TH2D hh = *h;  // materialize
        std::cout << "[INFO] truth TH2 entries = " << hh.GetEntries()
                  << ", integral(weighted) = " << hh.Integral() << "\n";
        DrawAndSave(hh, out_dir + "/pythia_truth_dr_vs_pair_pt.png",
                    "Pythia 8.3 truth single-b, pp #sqrt{s}=5.36 TeV (OS, weighted)");
    }

    // ======================= 2. pp24 DATA =======================
    {
        std::cout << "\n[INFO] pp24 data...\n";
        auto df = ROOT::RDataFrame("muon_pair_tree_sign2", std::vector<std::string>{pp24_file})
                      .Filter(data_cuts);
        ROOT::RDF::TH2DModel model("h_pp24", "pp24 data;p_{T}^{#mu#mu} [GeV];#DeltaR",
                                   nx, xbins.data(), ny, ylo, yhi);
        auto h = df.Histo2D(model, "pair_pt", "dr");
        TH2D hh = *h;
        std::cout << "[INFO] pp24 TH2 entries = " << hh.GetEntries()
                  << ", integral = " << hh.Integral() << "\n";
        DrawAndSave(hh, out_dir + "/pp24_data_dr_vs_pair_pt.png",
                    "pp 2024 data (2mu4), OS, raw counts");
    }

    // ======================= 3. PbPb DATA (23+24+25 combined) =======================
    {
        std::cout << "\n[INFO] PbPb data 2023+2024+2025...\n";
        auto df = ROOT::RDataFrame("muon_pair_tree_sign2", pbpb_files)
                      .Filter(data_cuts);
        ROOT::RDF::TH2DModel model("h_pbpb", "PbPb data 23+24+25;p_{T}^{#mu#mu} [GeV];#DeltaR",
                                   nx, xbins.data(), ny, ylo, yhi);
        auto h = df.Histo2D(model, "pair_pt", "dr");
        TH2D hh = *h;
        std::cout << "[INFO] PbPb TH2 entries = " << hh.GetEntries()
                  << ", integral = " << hh.Integral() << "\n";
        DrawAndSave(hh, out_dir + "/pbpb_data_dr_vs_pair_pt.png",
                    "PbPb 2023+2024+2025 data (single mu4), OS, raw counts");
    }

    std::cout << "\n[INFO] Done.\n";
}
