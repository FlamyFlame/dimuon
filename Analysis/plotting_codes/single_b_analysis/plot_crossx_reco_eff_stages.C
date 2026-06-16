// =============================================================================
// plot_crossx_reco_eff_stages.C
//
// Before/after efficiency-correction comparison for the single-b crossx, showing
// the impact of the reconstruction-efficiency PLACEHOLDER. Per pair-eta bin, the
// pair-pT differential cross-section is drawn as THREE lines:
//   - Raw (no correction)                         = *_corr_raw
//   - Reco-eff corrected (placeholder)            = *_corr_unfolded_reco
//   - Reco + trigger-eff corrected                = *_corr_unfolded_reco_trig
// (the "unfolded" stage is an identity placeholder, so raw == unfolded for now.)
//
// Uses the correction-stage histograms added by FillHistogramsCrossx (PP+PbPb).
// PbPb is summed over all centrality bins and years (combined, per analysis
// convention). pp24 is a single sample. See
// docs/tracking/reco_eff_placeholder_run2.md.
//
// Run: root -l -b -q plot_crossx_reco_eff_stages.C
// =============================================================================
#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../../Utilities/PbPbSampledLumi.h"

void plot_crossx_reco_eff_stages() {
    gStyle->SetOptStat(0);

    const std::string data_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data";
    const std::string out_dir  = data_dir + "/plots/sanity_check_crossx";
    gSystem->mkdir(out_dir.c_str(), true);

    const CommonEffcyConfig cfg{};
    const auto& eta_bins = cfg.pair_eta_proj_ranges_coarse_incl_gap;

    // The three correction stages to overlay (suffix, legend, color, marker).
    struct Stage { std::string suffix, label; int color, marker; };
    const std::vector<Stage> stages = {
        {"_corr_raw",               "Raw (no correction)",          kBlack,   24},
        {"_corr_unfolded_reco",     "Reco eff corr (placeholder)",  kBlue+1,  25},
        {"_corr_unfolded_reco_trig","Reco #times trig eff corr",    kRed+1,   20},
    };

    struct SampleSpec {
        std::string label;
        // PbPb: pattern-match per-centrality hists & sum; pp: exact base name.
        std::string pbpb_name_contains;  // "" for pp
        std::string pp_base_name;        // "" for PbPb
        std::string y_title;
        std::vector<std::pair<int, std::string>> year_paths;
        bool is_pp = false;
    };

    std::vector<SampleSpec> samples = {
        {"PbPb combined",
         "h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_ctr", "",
         "d#sigma/dp_{T} [nb GeV^{-1}]",
         {{23, data_dir + "/pbpb_2023"}, {24, data_dir + "/pbpb_2024"}, {25, data_dir + "/pbpb_2025"}}, false},
        {"PP 2024",
         "", "h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts",
         "d#sigma/dp_{T} [pb GeV^{-1}]",
         {{24, data_dir + "/pp_2024"}}, true},
    };

    for (const auto& spec : samples) {
        std::vector<std::pair<int,TFile*>> files; // (year, file) for luminosity-weighted combine
        for (const auto& [yr, dir] : spec.year_paths) {
            std::vector<std::string> candidates;
            if (spec.is_pp) {
                std::string base = dir + "/histograms_real_pairs_pp_20" + std::to_string(yr);
                candidates = { base + "_2mu4_nominal.root" };  // pp24 2mu4 crossx output
            } else {
                std::string base = dir + "/histograms_real_pairs_pbpb_20" + std::to_string(yr);
                candidates = { base + "_single_mu4_no_trg_plots_nominal.root" };  // PbPb single-mu4 crossx output
            }
            for (const auto& c : candidates) {
                TFile* f = TFile::Open(c.c_str(), "READ");
                if (f && !f->IsZombie()) { files.push_back({yr, f}); std::cout << "[INFO] Opened: " << c << "\n"; break; }
                if (f) { f->Close(); delete f; }
            }
        }
        if (files.empty()) { std::cerr << "[WARN] No files for " << spec.label << "\n"; continue; }

        // Combine a stage's TH2D across years. pp: single sample (simple). PbPb:
        // within each year sum the per-centrality stage histos, then combine years
        // by LUMINOSITY-WEIGHTED AVERAGE Sum(L_y·yr_sum)/Sum(L_y) (HF R_AA note Eq.3),
        // consistent with the crossx combined plotter and R_AA.
        auto getStage2D = [&](const Stage& st) -> TH2D* {
            TH2D* combined = nullptr;
            double sumL = 0.;
            for (auto& [yr, f] : files) {
                if (spec.is_pp) {
                    const std::string nm = spec.pp_base_name + st.suffix;
                    TH2D* h = dynamic_cast<TH2D*>(f->Get(nm.c_str()));
                    if (!h) continue;
                    if (!combined) { combined = dynamic_cast<TH2D*>(h->Clone((nm+"_comb").c_str())); combined->SetDirectory(nullptr); }
                    else combined->Add(h);
                } else {
                    // sum this year's per-centrality stage histos
                    TH2D* yr_sum = nullptr;
                    TIter next(f->GetListOfKeys()); TKey* key;
                    while ((key = (TKey*)next())) {
                        std::string nm = key->GetName();
                        const bool ends_with = nm.size() >= st.suffix.size() &&
                            nm.compare(nm.size()-st.suffix.size(), st.suffix.size(), st.suffix) == 0;
                        if (nm.find(spec.pbpb_name_contains) == std::string::npos || !ends_with) continue;
                        TH2D* h = dynamic_cast<TH2D*>(key->ReadObj());
                        if (!h) continue;
                        if (!yr_sum) { yr_sum = dynamic_cast<TH2D*>(h->Clone("pbpb_yr_sum")); yr_sum->SetDirectory(nullptr); }
                        else yr_sum->Add(h);
                    }
                    if (!yr_sum) continue;
                    const double L = PbPbMu4SampledLumiNb(yr);
                    if (!combined) { combined = dynamic_cast<TH2D*>(yr_sum->Clone("pbpb_stage_comb")); combined->SetDirectory(nullptr); combined->Scale(L); }
                    else combined->Add(yr_sum, L);
                    sumL += L;
                    delete yr_sum;
                }
            }
            if (!spec.is_pp && combined && sumL > 0.) combined->Scale(1.0 / sumL);
            return combined;
        };

        std::vector<TH2D*> h2(stages.size(), nullptr);
        bool ok = true;
        for (size_t s = 0; s < stages.size(); ++s) {
            h2[s] = getStage2D(stages[s]);
            if (!h2[s]) { std::cerr << "[WARN] Missing stage hist " << stages[s].suffix << " for " << spec.label << "\n"; ok = false; }
        }
        if (!ok) { for (auto& [yr,f] : files){f->Close();delete f;} continue; }

        int nrow = 2, ncol = 5;
        if ((int)eta_bins.size() <= 6) { nrow = 2; ncol = 3; }
        if ((int)eta_bins.size() <= 4) { nrow = 2; ncol = 2; }

        TCanvas c("c_reco_stages", (spec.label + " reco-eff stages").c_str(), 400*ncol, 350*nrow);
        c.Divide(ncol, nrow);
        std::vector<TH1D*> trash;

        for (size_t ieta = 0; ieta < eta_bins.size(); ++ieta) {
            c.cd((int)ieta + 1);
            gPad->SetLogx(); gPad->SetLogy();
            gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.13);

            const auto& eb = eta_bins[ieta];
            int y1 = h2[0]->GetYaxis()->FindBin(eb.first  + 1e-6);
            int y2 = h2[0]->GetYaxis()->FindBin(eb.second - 1e-6);

            std::vector<TH1D*> hp(stages.size(), nullptr);
            double ymax = 0;
            for (size_t s = 0; s < stages.size(); ++s) {
                hp[s] = h2[s]->ProjectionX(Form("st%zu_%zu_%d", s, ieta, rand()), y1, y2, "e");
                hp[s]->SetDirectory(nullptr);
                hp[s]->Scale(1.0, "width");
                hp[s]->SetLineColor(stages[s].color);
                hp[s]->SetMarkerColor(stages[s].color);
                hp[s]->SetMarkerStyle(stages[s].marker);
                hp[s]->SetMarkerSize(0.8);
                hp[s]->SetLineWidth(s == 0 ? 1 : 2);
                ymax = std::max(ymax, hp[s]->GetMaximum());
                trash.push_back(hp[s]);
            }
            hp[0]->SetMaximum(ymax * 3.0);
            hp[0]->SetMinimum(ymax * 1e-5);
            hp[0]->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
            hp[0]->GetYaxis()->SetTitle(spec.y_title.c_str());
            hp[0]->GetXaxis()->SetTitleSize(0.06); hp[0]->GetYaxis()->SetTitleSize(0.06);
            hp[0]->GetXaxis()->SetLabelSize(0.05); hp[0]->GetYaxis()->SetLabelSize(0.05);
            hp[0]->GetYaxis()->SetTitleOffset(1.45);
            hp[0]->SetTitle("");
            hp[0]->Draw("E1");
            for (size_t s = 1; s < stages.size(); ++s) hp[s]->Draw("E1 same");

            TLegend* leg = new TLegend(0.40, 0.66, 0.95, 0.92);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
            leg->AddEntry((TObject*)0, Form("#eta^{pair} #in [%.1f, %.1f]", eb.first, eb.second), "");
            for (size_t s = 0; s < stages.size(); ++s) leg->AddEntry(hp[s], stages[s].label.c_str(), "lpe");
            leg->Draw();
        }

        std::string safe = spec.label; std::replace(safe.begin(), safe.end(), ' ', '_');
        std::string out_path = out_dir + "/" + safe + "_reco_eff_stages_pair_pt_in_eta.png";
        c.SaveAs(out_path.c_str());
        std::cout << "[INFO] Saved: " << out_path << "\n";

        for (TH1D* h : trash) delete h;
        for (TH2D* h : h2) delete h;
        for (auto& [yr,f] : files) { f->Close(); delete f; }
    }
}
