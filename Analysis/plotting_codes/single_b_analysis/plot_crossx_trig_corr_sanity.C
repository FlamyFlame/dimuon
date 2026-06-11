#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"

void plot_crossx_trig_corr_sanity() {
    gStyle->SetOptStat(0);

    const std::string data_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data";
    const std::string out_dir  = data_dir + "/plots/sanity_check_crossx";
    gSystem->mkdir(out_dir.c_str(), true);

    const CommonEffcyConfig cfg{};
    const auto& eta_bins = cfg.pair_eta_proj_ranges_coarse_incl_gap;

    struct SampleSpec {
        std::string label;
        std::string h2_corrected;
        std::string h2_no_trig_corr;
        std::vector<std::pair<int, std::string>> year_paths;
        bool is_pp = false;
    };

    std::vector<SampleSpec> samples = {
        {"PbPb 10-20%",
         "h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_ctr10_20",
         "h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_ctr10_20_no_trig_corr",
         {{23, data_dir + "/pbpb_2023"}, {24, data_dir + "/pbpb_2024"}, {25, data_dir + "/pbpb_2025"}}
        },
        {"PP 2024",
         "h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts",
         "h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts_no_trig_corr",
         {{24, data_dir + "/pp_2024"}},
         true
        }
    };

    for (const auto& spec : samples) {
        std::vector<TFile*> files;
        for (const auto& [yr, dir] : spec.year_paths) {
            std::vector<std::string> candidates;
            if (spec.is_pp) {
                std::string base = dir + "/histograms_real_pairs_pp_20" + std::to_string(yr);
                const std::string trig = DatasetTriggerMap::Get(yr, 3);
                candidates = {
                    base + "_" + trig + "_nominal.root",
                    base + "_" + trig + "_coarse_q_eta_bin.root",
                };
            } else {
                std::string base = dir + "/histograms_real_pairs_pbpb_20" + std::to_string(yr);
                const std::string trig = DatasetTriggerMap::Get(yr, 1);
                candidates = {
                    base + "_" + trig + "_no_trg_plots_nominal.root",
                    base + "_" + trig + "_no_trg_plots_coarse_q_eta_bin.root",
                    base + "_" + trig + "_no_trg_plots_fine_q_eta_bin.root",
                };
            }
            for (const auto& c : candidates) {
                TFile* f = TFile::Open(c.c_str(), "READ");
                if (f && !f->IsZombie()) {
                    files.push_back(f);
                    std::cout << "[INFO] Opened: " << c << std::endl;
                    break;
                }
                if (f) { f->Close(); delete f; }
            }
        }
        if (files.empty()) {
            std::cerr << "[WARN] No files found for " << spec.label << std::endl;
            continue;
        }

        auto getCombined = [&](const std::string& name) -> TH2D* {
            TH2D* combined = nullptr;
            for (TFile* f : files) {
                TH2D* h = dynamic_cast<TH2D*>(f->Get(name.c_str()));
                if (!h) continue;
                if (!combined) {
                    combined = dynamic_cast<TH2D*>(h->Clone((name + "_comb").c_str()));
                    combined->SetDirectory(nullptr);
                } else {
                    combined->Add(h);
                }
            }
            return combined;
        };

        TH2D* h2_corr = getCombined(spec.h2_corrected);
        TH2D* h2_raw  = getCombined(spec.h2_no_trig_corr);

        if (!h2_corr) { std::cerr << "[WARN] Missing corrected hist: " << spec.h2_corrected << std::endl; continue; }
        if (!h2_raw)  { std::cerr << "[WARN] Missing no_trig_corr hist: " << spec.h2_no_trig_corr << std::endl; continue; }

        int nrow = 2, ncol = 5;
        if ((int)eta_bins.size() <= 6) { nrow = 2; ncol = 3; }
        if ((int)eta_bins.size() <= 4) { nrow = 2; ncol = 2; }

        TCanvas c("c_sanity", (spec.label + " raw vs trig-corr").c_str(), 400 * ncol, 350 * nrow);
        c.Divide(ncol, nrow);

        std::vector<TH1D*> trash;

        for (size_t ieta = 0; ieta < eta_bins.size(); ++ieta) {
            c.cd((int)ieta + 1);
            gPad->SetLogx();
            gPad->SetLogy();
            gPad->SetLeftMargin(0.16);
            gPad->SetBottomMargin(0.13);

            const auto& eb = eta_bins[ieta];
            int y1 = h2_corr->GetYaxis()->FindBin(eb.first  + 1e-6);
            int y2 = h2_corr->GetYaxis()->FindBin(eb.second - 1e-6);

            TH1D* hp_corr = h2_corr->ProjectionX(Form("corr_%zu_%d", ieta, rand()), y1, y2, "e");
            hp_corr->SetDirectory(nullptr);
            hp_corr->Scale(1.0, "width");

            TH1D* hp_raw = h2_raw->ProjectionX(Form("raw_%zu_%d", ieta, rand()), y1, y2, "e");
            hp_raw->SetDirectory(nullptr);
            hp_raw->Scale(1.0, "width");

            hp_raw->SetLineColor(kBlack);
            hp_raw->SetMarkerColor(kBlack);
            hp_raw->SetMarkerStyle(24);
            hp_raw->SetMarkerSize(0.8);
            hp_raw->SetLineWidth(1);

            hp_corr->SetLineColor(kRed + 1);
            hp_corr->SetMarkerColor(kRed + 1);
            hp_corr->SetMarkerStyle(20);
            hp_corr->SetMarkerSize(0.8);
            hp_corr->SetLineWidth(2);

            double ymax = std::max(hp_corr->GetMaximum(), hp_raw->GetMaximum());
            hp_corr->SetMaximum(ymax * 3.0);
            hp_corr->SetMinimum(ymax * 1e-5);
            hp_corr->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
            hp_corr->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb GeV^{-1}]");
            hp_corr->GetXaxis()->SetTitleSize(0.06);
            hp_corr->GetYaxis()->SetTitleSize(0.06);
            hp_corr->GetXaxis()->SetLabelSize(0.05);
            hp_corr->GetYaxis()->SetLabelSize(0.05);
            hp_corr->GetYaxis()->SetTitleOffset(1.45);
            hp_corr->SetTitle("");
            hp_corr->Draw("E1");
            hp_raw->Draw("E1 same");

            TLegend* leg = new TLegend(0.42, 0.68, 0.93, 0.92);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.042);
            leg->AddEntry((TObject*)0, Form("#eta^{pair} #in [%.1f, %.1f]", eb.first, eb.second), "");
            leg->AddEntry(hp_raw,  "Raw (no trig corr)", "lpe");
            leg->AddEntry(hp_corr, "No-corr trig eff", "lpe");
            leg->Draw();

            trash.push_back(hp_corr);
            trash.push_back(hp_raw);
        }

        std::string safe_label = spec.label;
        std::replace(safe_label.begin(), safe_label.end(), ' ', '_');
        std::string out_path = out_dir + "/" + safe_label + "_pair_pt_in_eta_subplots.png";
        c.SaveAs(out_path.c_str());
        std::cout << "[INFO] Saved: " << out_path << std::endl;

        for (TH1D* h : trash) delete h;
        delete h2_corr;
        delete h2_raw;
        for (TFile* f : files) { f->Close(); delete f; }
    }
}
