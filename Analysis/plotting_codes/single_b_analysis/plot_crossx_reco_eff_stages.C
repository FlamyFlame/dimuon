// =============================================================================
// plot_crossx_reco_eff_stages.C
//
// Correction-stage comparison for the single-b crossx. Per pair-eta bin, the pair-pT
// differential cross-section is drawn as THREE lines, TRIGGER FIRST (user, 2026-09-17;
// docs/tracking/pp24_trig_eff_hybrid_application.md §3f), so the trigger correction's own
// impact is the first step shown:
//   - Uncorrected                                       = *_corr_raw
//   - Trigger-efficiency corrected                      = *_corr_unfolded_trig
//   - Trigger- and reconstruction-efficiency corrected  = *_corr_unfolded_reco_trig
// FALLBACK, per sample: a histogram file produced BEFORE the trigger-first stage existed
// (2026-09-17, CorrectionStages.h) has no *_corr_unfolded_trig; such a sample is drawn with the
// older reco-first triple (raw / +reco / +reco+trig) and the fallback is PRINTED. The final
// stage is the same either way (the corrections commute). Pb+Pb on disk is in that state until
// its next crossx refill.
// The reconstruction efficiency is the measured pp24-fullsim 3D PAIR efficiency for pp and the
// Run-2 single-muon PLACEHOLDER for PbPb -- the legend says so for Pb+Pb.
// (the "unfolded" stage is an identity placeholder, so raw == unfolded for now.)
//
// Uses the correction-stage histograms added by FillHistogramsCrossx (PP+PbPb).
// PbPb is summed over all centrality bins and years (combined, per analysis
// convention). pp24 is a single sample. See
// docs/tracking/reco_eff_placeholder_run2.md.
//
// Muon working point: `use_tight_wp` (default TRUE = the nominal Tight, unsuffixed histogram
// files; false = the Medium WP-systematic `_medium_wp` files, docs/muon_wp_registry.md). The
// figure states the dataset, trigger and working point in its legend strip.
//
// Run: root -l -b -q plot_crossx_reco_eff_stages.C                 (both samples, Tight)
//      root -l -b -q 'plot_crossx_reco_eff_stages.C(false)'        (pp only -- the pp pipeline's Stage 7)
//      root -l -b -q 'plot_crossx_reco_eff_stages.C(false,false)'  (pp only, Medium WP)
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
#include <cmath>
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../../Utilities/CommonLogYRange.h"
#include "../../Utilities/PbPbSampledLumi.h"
#include "../../Utilities/PairEtaPanelBins.h"

void plot_crossx_reco_eff_stages(bool include_pbpb = true, bool use_tight_wp = true) {
    gStyle->SetOptStat(0);
    const std::string wp_suffix = use_tight_wp ? "" : "_medium_wp";
    const std::string wp_text   = use_tight_wp ? "tight WP" : "medium WP";

    const std::string data_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data";
    const std::string out_dir  = data_dir + "/plots/sanity_check_crossx";
    gSystem->mkdir(out_dir.c_str(), true);

    const CommonEffcyConfig cfg{};
    const auto& eta_bins = cfg.pair_eta_proj_ranges_coarse_incl_gap;

    // The three correction stages to overlay (suffix, legend, color, marker), TRIGGER FIRST.
    // The reconstruction-efficiency label is SAMPLE-DEPENDENT. pp applies the measured
    // pp24-fullsim 3D pair efficiency; Pb+Pb still applies the Run 2 single-muon PLACEHOLDER, and
    // a reader of the Pb+Pb figure must be able to see that from the figure itself.
    struct Stage { std::string suffix, label, label_pbpb; int color, marker; };
    const std::vector<Stage> stages_trig_first = {
        {"_corr_raw",               "Uncorrected", "Uncorrected", kBlack, 24},
        {"_corr_unfolded_trig",     "+ trigger efficiency", "+ trigger efficiency", kBlue+1, 25},
        {"_corr_unfolded_reco_trig","+ reconstruction efficiency",
                                    "+ reconstruction efficiency (Run 2 placeholder)", kRed+1, 20},
    };
    // The pre-2026-09-17 reco-first triple, used ONLY for a file that lacks the trigger-first
    // stage (see the header). Same first and last stage.
    const std::vector<Stage> stages_reco_first = {
        {"_corr_raw",               "Uncorrected", "Uncorrected", kBlack, 24},
        {"_corr_unfolded_reco",     "+ reconstruction efficiency",
                                    "+ reconstruction efficiency (Run 2 placeholder)", kBlue+1, 25},
        {"_corr_unfolded_reco_trig","+ trigger efficiency", "+ trigger efficiency", kRed+1, 20},
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
         {{23, data_dir + "/pbpb_2023"}, {24, data_dir + "/pbpb_2024"},
          {25, data_dir + "/pbpb_2025"}, {26, data_dir + "/pbpb_2026"}}, false},
        {"PP 2024",
         "", "h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts",
         "d#sigma/dp_{T} [pb GeV^{-1}]",
         {{24, data_dir + "/pp_2024"}}, true},
    };

    for (const auto& spec : samples) {
        if (!include_pbpb && !spec.is_pp) {
            std::cout << "SKIP (include_pbpb=false): " << spec.label << std::endl;
            continue;
        }
        std::vector<TH2D*> h2;
        std::vector<Stage> stages;   // the triple actually drawn for THIS sample
        std::vector<std::pair<int,TFile*>> files; // (year, file) for luminosity-weighted combine
        for (const auto& [yr, dir] : spec.year_paths) {
            std::vector<std::string> candidates;
            if (spec.is_pp) {
                std::string base = dir + "/histograms_real_pairs_pp_20" + std::to_string(yr);
                candidates = { base + "_2mu4_nominal" + wp_suffix + ".root" };  // pp24 2mu4 crossx output
            } else {
                std::string base = dir + "/histograms_real_pairs_pbpb_20" + std::to_string(yr);
                candidates = { base + "_single_mu4_no_trg_plots_nominal" + wp_suffix + ".root" };  // PbPb single-mu4 crossx output
            }
            for (const auto& c : candidates) {
                // Probe before opening: TFile::Open on a non-existent path prints a raw
                // "Error in <TFile::TFile>: file ... does not exist" to stderr, which for a
                // year that simply has not been produced yet is noise, not an error.
                if (gSystem->AccessPathName(c.c_str())) continue;
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

        // Trigger-first only if EVERY opened file carries that stage; otherwise the older
        // reco-first triple for the whole sample. Decided per FILE, not on the combined
        // histogram: getStage2D silently skips a year that lacks the stage, so a mixed set of
        // years (some refilled after 2026-09-17, some not) would otherwise draw an intermediate
        // curve summed over a SUBSET of the years against first/last stages summed over all of
        // them -- a year-set mismatch that would read as a correction step in the ratio pad.
        {
            std::vector<int> missing;
            for (auto& [yr, f] : files) {
                bool has = false;
                if (spec.is_pp) {
                    has = f->Get((spec.pp_base_name + stages_trig_first[1].suffix).c_str()) != nullptr;
                } else {
                    TIter next(f->GetListOfKeys()); TKey* key;
                    while ((key = (TKey*)next())) {
                        const std::string nm = key->GetName();
                        const std::string& suf = stages_trig_first[1].suffix;
                        if (nm.find(spec.pbpb_name_contains) != std::string::npos &&
                            nm.size() >= suf.size() &&
                            nm.compare(nm.size() - suf.size(), suf.size(), suf) == 0) { has = true; break; }
                    }
                }
                if (!has) missing.push_back(yr);
            }
            if (missing.empty()) stages = stages_trig_first;
            else {
                stages = stages_reco_first;
                std::cout << "[INFO] " << spec.label << ": no '" << stages_trig_first[1].suffix
                          << "' stage in the histogram file(s) of year(s)";
                for (int yr : missing) std::cout << " " << yr;
                std::cout << " (produced before 2026-09-17) -> drawing the reco-first triple "
                             "raw / +reco / +reco+trig for the whole sample instead" << std::endl;
            }
        }
        h2.assign(stages.size(), nullptr);
        bool ok = true;
        for (size_t s = 0; s < stages.size(); ++s) {
            h2[s] = getStage2D(stages[s]);
            if (!h2[s]) { std::cerr << "[WARN] Missing stage hist " << stages[s].suffix << " for " << spec.label << "\n"; ok = false; }
        }
        if (!ok) { for (auto& [yr,f] : files){f->Close();delete f;} continue; }

        // nrows >= ncols, nrows ~ sqrt(N)  (subplot-layout convention). Nine pair-eta bins are
        // 3x3, never 5x2 -- the sibling trigger-corrections figure uses the same grid, and two
        // figures of the same spectrum must not be laid out differently.
        const int n_eta = (int)eta_bins.size();
        const int ncol = (n_eta <= 3) ? n_eta : (int)std::ceil(std::sqrt((double)n_eta));
        const int nrow = (int)std::ceil((double)n_eta / ncol);

        TCanvas c("c_reco_stages", (spec.label + " reco-eff stages").c_str(), 430*ncol, 380*nrow);
        c.Divide(ncol, nrow);
        std::vector<TH1D*> trash;
        std::vector<std::vector<TH1D*>> panels(eta_bins.size());

        // PASS 1 — build every panel's stage curves WITHOUT drawing, so the common
        // log-y range can be derived from all of them (Utilities/CommonLogYRange.h).
        for (size_t ieta = 0; ieta < eta_bins.size(); ++ieta) {
            const auto& eb = eta_bins[ieta];
            const auto yb = PairEtaPanels::Bins(h2[0]->GetYaxis(), eb,
                                               "plot_crossx_reco_eff_stages");
            int y1 = yb.first, y2 = yb.second;

            std::vector<TH1D*> hp(stages.size(), nullptr);
            for (size_t s = 0; s < stages.size(); ++s) {
                hp[s] = h2[s]->ProjectionX(Form("st%zu_%zu_%d", s, ieta, rand()), y1, y2, "e");
                hp[s]->SetDirectory(nullptr);
                hp[s]->Scale(1.0, "width");
                hp[s]->SetLineColor(stages[s].color);
                hp[s]->SetMarkerColor(stages[s].color);
                hp[s]->SetMarkerStyle(stages[s].marker);
                hp[s]->SetMarkerSize(0.8);
                hp[s]->SetLineWidth(s == 0 ? 1 : 2);
                trash.push_back(hp[s]);
            }
            hp[0]->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
            hp[0]->GetYaxis()->SetTitle(spec.y_title.c_str());
            hp[0]->GetXaxis()->SetTitleSize(0.06); hp[0]->GetYaxis()->SetTitleSize(0.06);
            hp[0]->GetXaxis()->SetLabelSize(0.05); hp[0]->GetYaxis()->SetLabelSize(0.05);
            hp[0]->GetYaxis()->SetTitleOffset(1.55);
            hp[0]->SetTitle("");
            panels[ieta] = hp;
        }

        // ONE log-y scale for the whole PNG, low enough to keep every non-empty
        // point of every panel and every correction stage inside the frame. The
        // previous fixed floor (ymax * 1e-5) cropped the highest-pair-pT points,
        // whose spread below the panel maximum exceeds 1e5.
        {
            std::vector<TH1*> flat;
            for (const auto& hp : panels)
                for (TH1D* h : hp) if (h) flat.push_back(h);
            ApplyCommonLogYRange(flat, 3.0);  // a little headroom above the highest point
        }

        // PASS 2 — draw. Each panel is split into a spectrum pad and a ratio pad, and the
        // legend gets its OWN strip at the top of the canvas rather than a box inside a frame:
        // widening it to fit the longest entry only traded clipping for an overlap with the
        // highest-pair-pT points (atlas-plotting.md -- reserve space, do not overlay).
        const double head = 0.085;
        TPad* head_pad = new TPad("head_pad", "", 0.0, 1.0 - head, 1.0, 1.0);
        TPad* body_pad = new TPad("body_pad", "", 0.0, 0.0, 1.0, 1.0 - head);
        head_pad->SetFillStyle(0); body_pad->SetFillStyle(0);
        head_pad->Draw(); body_pad->Draw();

        head_pad->cd();
        {
            TLegend* leg = new TLegend(0.05, 0.05, 0.98, 0.95);
            leg->SetNColumns((int)stages.size());
            leg->SetBorderSize(0); leg->SetFillStyle(0);
            // The Pb+Pb reconstruction entry carries "(Run 2 placeholder)" and is much longer
            // than the pp one; at a single shared text size its closing bracket ran into the
            // next column's marker.
            leg->SetTextSize(spec.is_pp ? 0.32 : 0.24);
            for (size_t s = 0; s < stages.size(); ++s)
                leg->AddEntry(panels[0][s],
                              (spec.is_pp ? stages[s].label : stages[s].label_pbpb).c_str(), "lpe");
            leg->Draw();
        }

        // ONE ratio range for all nine panels. With a per-panel range the same visual height
        // means a different correction from one panel to the next, which is exactly the
        // comparison the pad exists to make. Built first, over every panel's ratios.
        std::vector<std::vector<TH1D*>> all_ratios(eta_bins.size());
        double rmin_all = 1.0, rmax_all = 1.0;
        for (size_t ieta = 0; ieta < eta_bins.size(); ++ieta) {
            auto& hp = panels[ieta];
            if (hp.empty() || !hp[0]) continue;
            for (size_t st = 1; st < stages.size(); ++st) {
                TH1D* r = (TH1D*)hp[st]->Clone(Form("r%zu_%zu", st, ieta));
                r->SetDirectory(nullptr);
                r->Divide(hp[0]);
                trash.push_back(r);
                all_ratios[ieta].push_back(r);
                for (int b = 1; b <= r->GetNbinsX(); ++b) {
                    if (hp[0]->GetBinContent(b) <= 0. || r->GetBinContent(b) <= 0.) continue;
                    rmin_all = std::min(rmin_all, r->GetBinContent(b));
                    rmax_all = std::max(rmax_all, r->GetBinContent(b));
                }
            }
        }
        std::cout << "[INFO] shared ratio range for " << spec.label << ": ["
                  << 0.9 * rmin_all << ", " << 1.1 * rmax_all << "]" << std::endl;

        body_pad->cd();
        body_pad->Divide(ncol, nrow);
        for (size_t ieta = 0; ieta < eta_bins.size(); ++ieta) {
            auto& hp = panels[ieta];
            if (hp.empty() || !hp[0]) continue;

            body_pad->cd((int)ieta + 1);
            const auto& eb = eta_bins[ieta];

            // Upper pad: the three spectra. Lower pad: corrected / uncorrected -- the size of
            // the correction is the point of the figure and cannot be read off a 7-decade log
            // axis by eye (R3).
            TPad* up = new TPad(Form("up%zu", ieta), "", 0, 0.38, 1, 1);
            TPad* lo = new TPad(Form("lo%zu", ieta), "", 0, 0.0,  1, 0.38);
            up->SetLogx(); up->SetLogy();
            up->SetLeftMargin(0.19); up->SetBottomMargin(0.02); up->SetTopMargin(0.04);
            lo->SetLogx(); lo->SetGridy();
            lo->SetLeftMargin(0.19); lo->SetTopMargin(0.02); lo->SetBottomMargin(0.30);
            up->Draw(); lo->Draw();

            up->cd();
            hp[0]->GetXaxis()->SetLabelSize(0);
            hp[0]->GetXaxis()->SetTitleSize(0);
            hp[0]->GetYaxis()->SetTitleSize(0.075);
            hp[0]->GetYaxis()->SetLabelSize(0.065);
            hp[0]->GetYaxis()->SetTitleOffset(1.20);
            hp[0]->Draw("E1");
            for (size_t s = 1; s < stages.size(); ++s) hp[s]->Draw("E1 same");
            TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.075);
            t.DrawLatex(0.24, 0.10, Form("#eta^{pair} #in [%.1f, %.1f]", eb.first, eb.second));
            // Dataset, trigger and working point, top-right of every spectrum pad (empty on a
            // steeply falling log-y spectrum), as the crossx panels state them.
            t.SetTextSize(0.058); t.SetTextAlign(33);
            t.DrawLatex(0.89, 0.91, (spec.is_pp ? "pp 2024, 2mu4, " + wp_text
                                                : "Pb+Pb combined, mu4, " + wp_text).c_str());

            lo->cd();
            std::vector<TH1D*>& ratios = all_ratios[ieta];
            for (size_t i = 0; i < ratios.size(); ++i) {
                TH1D* r = ratios[i];
                r->SetTitle("");
                r->GetYaxis()->SetTitle("corrected / uncorrected");
                // The reference line at 1 must be inside the frame (a 1/eps correction is >= 1,
                // so the lower edge is pulled down to just below 1).
                r->GetYaxis()->SetRangeUser(std::min(0.9 * rmin_all, 0.9), 1.1 * rmax_all);
                r->GetYaxis()->SetNdivisions(505);
                r->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
                r->GetXaxis()->SetTitleSize(0.120); r->GetXaxis()->SetLabelSize(0.105);
                r->GetYaxis()->SetTitleSize(0.090); r->GetYaxis()->SetLabelSize(0.100);
                r->GetYaxis()->SetTitleOffset(0.95);
                r->GetXaxis()->SetTitleOffset(1.00);
                r->Draw(i == 0 ? "E1" : "E1 same");
            }
            if (!ratios.empty()) {
                // The "no correction" reference. R3: a ratio pad needs a line at 1.
                TLine* one = new TLine(ratios[0]->GetXaxis()->GetXmin(), 1.0,
                                       ratios[0]->GetXaxis()->GetXmax(), 1.0);
                one->SetLineStyle(2);
                one->SetLineColor(kGray + 2);
                one->Draw("SAME");
            }
        }

        std::string safe = spec.label; std::replace(safe.begin(), safe.end(), ' ', '_');
        std::string out_path = out_dir + "/" + safe + "_reco_eff_stages_pair_pt_in_eta" + wp_suffix + ".png";
        c.SaveAs(out_path.c_str());
        std::cout << "[INFO] Saved: " << out_path << "\n";

        for (TH1D* h : trash) delete h;
        for (TH2D* h : h2) delete h;
        for (auto& [yr,f] : files) { f->Close(); delete f; }
    }
}
