// H1 investigation: Compare fitted TF1 vs raw binned efficiency from Pipeline 2
// Standalone macro — does NOT modify any analysis code.
// Saves plots to investigation/ subdirectory.
// To undo: simply delete this macro and the investigation/ directory.

#include "TFile.h"
#include "TF1.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TKey.h"
#include <vector>
#include <string>
#include <map>
#include <cmath>

void investigate_H1_fit_vs_data() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    std::string base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/";
    std::string fit_path = base + "trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit.root";
    std::string hist_path = base + "histograms_real_pairs_pbpb_2023_single_mu4_coarse_q_eta_bin.root";
    std::string outdir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pbpb_trigger_efficiency/mu4/investigation/H1_fit_vs_data/";
    gSystem->mkdir(outdir.c_str(), true);

    auto f_fit = TFile::Open(fit_path.c_str());
    auto f_hist = TFile::Open(hist_path.c_str());
    if (!f_fit || !f_hist) { std::cerr << "Cannot open files" << std::endl; return; }

    // Load all TF1s into a map
    std::map<std::string, TF1*> tf1_map;
    TIter nextk(f_fit->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextk())) {
        TObject* obj = key->ReadObj();
        if (obj->InheritsFrom("TF1")) {
            tf1_map[obj->GetName()] = (TF1*)obj;
        }
    }
    printf("Loaded %d TF1s\n", (int)tf1_map.size());

    // q_eta bins (fine, gap-excluded) — must match what EvaluateSingleMuonEffcyPtFitted uses
    struct QEtaBin { float lo, hi; std::string suffix; };
    std::vector<QEtaBin> qeta_bins = {
        {-2.4f, -2.0f, "minus2_40_TO_minus2_00"},
        {-2.0f, -1.6f, "minus2_00_TO_minus1_60"},
        {-1.6f, -1.3f, "minus1_60_TO_minus1_30"},
        {-0.9f, -0.5f, "minus0_90_TO_minus0_50"},
        {-0.5f, -0.1f, "minus0_50_TO_minus0_10"},
        {0.1f, 0.5f, "0_10_TO_0_50"},
        {0.5f, 1.0f, "0_50_TO_1_00"},
        {1.3f, 1.6f, "1_30_TO_1_60"},
        {1.6f, 2.0f, "1_60_TO_2_00"},
        {2.0f, 2.2f, "2_00_TO_2_20"}
    };

    struct CtrBin { std::string label; std::string suffix; };
    std::vector<CtrBin> ctr_bins = {
        {"0-5%", "_ctr0_5"}, {"5-10%", "_ctr5_10"}, {"10-20%", "_ctr10_20"},
        {"20-30%", "_ctr20_30"}, {"30-50%", "_ctr30_50"}, {"50-80%", "_ctr50_80"}
    };

    std::vector<std::pair<std::string, std::string>> musigns = {
        {"_sign1", "mu+"},
        {"_sign2", "mu-"}
    };

    // Summary statistics
    int n_bins_total = 0, n_bins_undershoots = 0, n_bins_overshoots = 0;
    double sum_mismatch = 0;
    std::map<std::string, double> ctr_avg_mismatch;
    std::map<std::string, int> ctr_n_bins;

    for (const auto& ctr : ctr_bins) {
        ctr_avg_mismatch[ctr.suffix] = 0;
        ctr_n_bins[ctr.suffix] = 0;

        for (const auto& [msign, mlabel] : musigns) {
            // Get the 2D divided histogram
            std::string h2d_name = "h_pt2nd_vs_q_eta2nd" + ctr.suffix + msign + "_2mu4_sepr_divided";
            auto h2d = (TH2D*)f_hist->Get(h2d_name.c_str());
            if (!h2d) { printf("MISSING: %s\n", h2d_name.c_str()); continue; }

            // Also get the denominator (total counts) to understand statistics
            std::string h2d_denom_name = "h_pt2nd_vs_q_eta2nd" + ctr.suffix + msign + "_mu4_sepr";
            auto h2d_denom = (TH2D*)f_hist->Get(h2d_denom_name.c_str());

            for (const auto& qeta : qeta_bins) {
                std::string tf1_name = "f_pt2nd_vs_q_eta2nd" + ctr.suffix + msign + "_2mu4_sepr_py_" + qeta.suffix + "_divided";
                auto it = tf1_map.find(tf1_name);
                if (it == tf1_map.end()) continue;
                TF1* fit = it->second;

                // Project 2D hist onto pT for this q_eta range
                int bin_lo = h2d->GetXaxis()->FindBin(qeta.lo + 0.001);
                int bin_hi = h2d->GetXaxis()->FindBin(qeta.hi - 0.001);
                TH1D* h_proj = h2d->ProjectionY(Form("proj_%s_%s_%s", ctr.suffix.c_str(), msign.c_str(), qeta.suffix.c_str()), bin_lo, bin_hi);

                // Also project denominator for error estimation
                TH1D* h_denom_proj = nullptr;
                if (h2d_denom) {
                    h_denom_proj = h2d_denom->ProjectionY(Form("denom_%s_%s_%s", ctr.suffix.c_str(), msign.c_str(), qeta.suffix.c_str()), bin_lo, bin_hi);
                }

                // Compare fit vs data for pT > 5 GeV (above threshold)
                for (int iy = 1; iy <= h_proj->GetNbinsX(); iy++) {
                    double pt = h_proj->GetBinCenter(iy);
                    if (pt < 5.0) continue;
                    double data_val = h_proj->GetBinContent(iy);
                    if (data_val < 0.01 || data_val > 1.0) continue;
                    double n_denom = h_denom_proj ? h_denom_proj->GetBinContent(iy) : 1;
                    if (n_denom < 5) continue; // skip bins with too few events

                    double fit_val = fit->Eval(pt);
                    if (fit_val < 0.01) fit_val = 0.01;
                    if (fit_val > 1.0) fit_val = 1.0;

                    double mismatch = (fit_val - data_val) / data_val;
                    n_bins_total++;
                    if (fit_val < data_val) n_bins_undershoots++;
                    else n_bins_overshoots++;
                    sum_mismatch += mismatch;
                    ctr_avg_mismatch[ctr.suffix] += mismatch;
                    ctr_n_bins[ctr.suffix]++;
                }

                // Make overlay plot for a representative selection
                if (ctr.suffix == "_ctr10_20" && msign == "_sign1") {
                    auto c = new TCanvas("c", "", 800, 600);
                    c->SetLeftMargin(0.12);
                    c->SetLogx();
                    h_proj->SetMarkerStyle(20);
                    h_proj->SetMarkerSize(0.7);
                    h_proj->GetXaxis()->SetTitle("p_{T} [GeV]");
                    h_proj->GetYaxis()->SetTitle("#varepsilon(2mu4 | mu4, #DeltaR > 0.8)");
                    h_proj->GetYaxis()->SetRangeUser(0, 1.1);
                    h_proj->GetXaxis()->SetRangeUser(3.5, 60);
                    h_proj->SetTitle("");
                    h_proj->Draw("E");
                    fit->SetLineColor(kRed);
                    fit->SetLineWidth(2);
                    fit->Draw("same");
                    auto leg = new TLegend(0.55, 0.2, 0.88, 0.4);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->AddEntry(h_proj, "Data (binned)", "lp");
                    leg->AddEntry(fit, "Fermi+log fit", "l");
                    leg->Draw();
                    TLatex tex;
                    tex.SetNDC();
                    tex.SetTextSize(0.03);
                    tex.DrawLatex(0.15, 0.88, Form("%s, %s, q#eta [%.1f, %.1f]", ctr.label.c_str(), mlabel.c_str(), qeta.lo, qeta.hi));
                    c->SaveAs(Form("%sfit_vs_data_%s_%s_%s.png", outdir.c_str(), ctr.suffix.c_str()+1, msign.c_str()+1, qeta.suffix.c_str()));
                    delete c;
                }

                delete h_proj;
                if (h_denom_proj) delete h_denom_proj;
            }
        }
    }

    // Print summary
    printf("\n=== H1 SUMMARY: Fit vs Data Mismatch ===\n");
    printf("Total pT bins compared (pT>5, N_denom>=5): %d\n", n_bins_total);
    printf("Bins where fit < data (undershoots): %d (%.1f%%)\n", n_bins_undershoots, 100.0*n_bins_undershoots/n_bins_total);
    printf("Bins where fit > data (overshoots):  %d (%.1f%%)\n", n_bins_overshoots, 100.0*n_bins_overshoots/n_bins_total);
    printf("Average relative mismatch (fit-data)/data: %.4f\n", sum_mismatch/n_bins_total);
    printf("\nPer centrality:\n");
    for (const auto& ctr : ctr_bins) {
        int n = ctr_n_bins[ctr.suffix];
        double avg = n > 0 ? ctr_avg_mismatch[ctr.suffix] / n : 0;
        printf("  %s: avg mismatch = %+.4f (%d bins)\n", ctr.label.c_str(), avg, n);
    }

    f_fit->Close();
    f_hist->Close();
    printf("\nPlots saved to %s\n", outdir.c_str());
}
