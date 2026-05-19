#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TSystem.h"
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <cmath>

bool isBNL = true;

std::string bnl_base() { return "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/"; }
std::string mac_base() { return "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/"; }
std::string base_dir() { return isBNL ? bnl_base() : mac_base(); }
std::string plot_base() { return isBNL ? bnl_base() + "plots/" : mac_base() + "plots/"; }

struct PlateauResult {
    double value;
    double error;
    int npoints;
};

PlateauResult extractPlateau(TGraphAsymmErrors* g, double dR_lo = 1.0, double dR_hi = 3.0) {
    double sum_w = 0, sum_wy = 0;
    int npts = 0;
    for (int i = 0; i < g->GetN(); i++) {
        double x, y;
        g->GetPoint(i, x, y);
        if (x < dR_lo || x > dR_hi || y <= 0.01) continue;
        double ey = std::max(g->GetErrorYhigh(i), g->GetErrorYlow(i));
        if (ey <= 0) ey = 0.01;
        double w = 1.0 / (ey * ey);
        sum_w += w;
        sum_wy += w * y;
        npts++;
    }
    if (npts == 0 || sum_w == 0) return {1.0, 0.0, 0};
    return {sum_wy / sum_w, 1.0 / std::sqrt(sum_w), npts};
}

TGraphAsymmErrors* normalizeGraph(TGraphAsymmErrors* g, double plateau) {
    auto gn = (TGraphAsymmErrors*)g->Clone(Form("%s_platnorm", g->GetName()));
    for (int i = 0; i < gn->GetN(); i++) {
        double x, y;
        gn->GetPoint(i, x, y);
        gn->SetPoint(i, x, y / plateau);
        gn->SetPointEYhigh(i, gn->GetErrorYhigh(i) / plateau);
        gn->SetPointEYlow(i, gn->GetErrorYlow(i) / plateau);
    }
    return gn;
}

void plateau_normalize_dR_corr() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    struct YearInfo { int year; std::string path; int color; int marker; };
    std::vector<YearInfo> years = {
        {23, base_dir() + "pbpb_2023/histograms_real_pairs_pbpb_2023_single_mu4_coarse_q_eta_bin.root", kBlue, 20},
        {24, base_dir() + "pbpb_2024/histograms_real_pairs_pbpb_2024_single_mu4_coarse_q_eta_bin.root", kRed, 21},
        {25, base_dir() + "pbpb_2025/histograms_real_pairs_pbpb_2025_single_mu4_coarse_q_eta_bin.root", kGreen+2, 22}
    };

    std::map<int, TFile*> files;
    for (auto& y : years) {
        files[y.year] = TFile::Open(y.path.c_str());
        if (!files[y.year] || files[y.year]->IsZombie()) {
            std::cerr << "Cannot open " << y.path << std::endl;
            return;
        }
    }

    std::vector<std::string> ctr_cats = {""};
    std::vector<std::string> ctr_labels = {""};
    std::vector<std::pair<std::string, std::string>> ctr_bins = {
        {"0-5%", "_ctr0_5"}, {"5-10%", "_ctr5_10"}, {"10-20%", "_ctr10_20"},
        {"20-30%", "_ctr20_30"}, {"30-50%", "_ctr30_50"}, {"50-80%", "_ctr50_80"}
    };
    for (auto& [lab, suf] : ctr_bins) {
        ctr_cats.push_back(suf);
        ctr_labels.push_back(lab);
    }

    struct PtSlice { std::string suffix; std::string label; };
    std::vector<PtSlice> pt_slices = {
        {"_pt8_12", "8 < p_{T}^{pair} < 12 GeV"},
        {"_pt12_20", "12 < p_{T}^{pair} < 20 GeV"},
        {"_pt20_40", "20 < p_{T}^{pair} < 40 GeV"},
        {"_pt40_120", "40 < p_{T}^{pair} < 120 GeV"}
    };

    struct PlotVar {
        std::string dr_var;
        std::string xtitle;
        double xmin, xmax;
    };
    std::vector<PlotVar> plot_vars = {
        {"DR",        "#DeltaR",  0.0, 5.75},
        {"DR_zoomin", "#DeltaR",  0.0, 0.8},
        {"DR_0_2",    "#DeltaR",  0.0, 0.2}
    };

    std::string pair_suffix = "_pair_mu4";

    // ============================================================
    // MODE 1: pT-integrated plateau normalization
    // ============================================================
    auto makeNormPlots = [&](const std::string& mode_label, const std::string& outdir_suffix,
                              bool use_pt_binned_plateau) {
        std::string base_outdir = plot_base() + "pbpb_trigger_efficiency/mu4/dR_pair_level_" + outdir_suffix + "/";

        // Also save normalized graphs to ROOT files
        std::map<int, TFile*> out_files;
        for (auto& y : years) {
            std::string root_outdir = base_dir() + "pbpb_20" + std::to_string(y.year) + "/";
            std::string outpath = root_outdir + "dR_corr_plateau_norm_" + outdir_suffix + "_pbpb_20" + std::to_string(y.year) + ".root";
            out_files[y.year] = TFile::Open(outpath.c_str(), "RECREATE");
        }

        for (const auto& sign : {"_op", "_ss"}) {
            std::string sign_label = (std::string(sign) == "_op") ? "OS" : "SS";

            for (int icat = 0; icat < (int)ctr_cats.size(); icat++) {
                std::string cat = ctr_cats[icat];
                std::string cat_label = ctr_labels[icat];
                std::string ctr_subdir = cat.empty() ? "ctr_integrated" : "ctr_binned";

                // Extract pT-integrated plateau from full-range DR graph (per year)
                std::map<int, double> pt_int_plateau;
                for (auto& y : years) {
                    std::string gname = "g_DR" + std::string(sign) + pair_suffix + "_invw_num" + cat + "_divided";
                    auto g = (TGraphAsymmErrors*)files[y.year]->Get(gname.c_str());
                    if (!g) { pt_int_plateau[y.year] = 1.0; continue; }
                    auto plat = extractPlateau(g, 1.0, 3.0);
                    pt_int_plateau[y.year] = plat.npoints > 0 ? plat.value : 1.0;
                    std::cout << Form("  %s yr%d %s%s: plateau=%.4f +/- %.4f (%d pts)",
                        mode_label.c_str(), y.year, sign, cat.c_str(), plat.value, plat.error, plat.npoints) << std::endl;
                }

                // Extract pT-binned plateaus (for mode 2)
                std::map<int, std::map<std::string, double>> pt_bin_plateaus;
                if (use_pt_binned_plateau) {
                    for (auto& y : years) {
                        for (auto& sl : pt_slices) {
                            std::string gname = "g_DR" + std::string(sign) + pair_suffix + "_invw_num" + cat + sl.suffix + "_divided";
                            auto g = (TGraphAsymmErrors*)files[y.year]->Get(gname.c_str());
                            if (!g) { pt_bin_plateaus[y.year][sl.suffix] = pt_int_plateau[y.year]; continue; }
                            auto plat = extractPlateau(g, 1.0, 3.0);
                            pt_bin_plateaus[y.year][sl.suffix] = plat.npoints >= 3 ? plat.value : pt_int_plateau[y.year];
                            if (plat.npoints >= 3) {
                                std::cout << Form("    %s yr%d %s%s%s: plateau=%.4f (%d pts)",
                                    mode_label.c_str(), y.year, sign, cat.c_str(), sl.suffix.c_str(), plat.value, plat.npoints) << std::endl;
                            }
                        }
                    }
                }

                // --- Plot pT-integrated dR corrections (normalized) ---
                for (const auto& pv : plot_vars) {
                    std::string gname_base = "g_" + pv.dr_var + std::string(sign) + pair_suffix + "_invw_num" + cat + "_divided";

                    auto c = new TCanvas("c", "", 800, 600);
                    c->SetLeftMargin(0.12);
                    c->SetRightMargin(0.05);
                    auto frame = c->DrawFrame(pv.xmin, 0.85, pv.xmax, 1.10);
                    frame->GetXaxis()->SetTitle(pv.xtitle.c_str());
                    frame->GetYaxis()->SetTitle("#varepsilon_{#DeltaR}^{pair} (plateau normalized)");

                    auto leg = new TLegend(0.55, 0.72, 0.90, 0.90);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);

                    bool any_drawn = false;
                    for (auto& y : years) {
                        auto g = (TGraphAsymmErrors*)files[y.year]->Get(gname_base.c_str());
                        if (!g) continue;
                        double plat = pt_int_plateau[y.year];
                        auto gn = normalizeGraph(g, plat);
                        gn->SetMarkerColor(y.color);
                        gn->SetLineColor(y.color);
                        gn->SetMarkerStyle(y.marker);
                        gn->SetMarkerSize(0.8);
                        gn->Draw("P same");
                        leg->AddEntry(gn, Form("PbPb 20%d (plat=%.3f)", y.year, plat), "lp");

                        out_files[y.year]->cd();
                        std::string save_name = "g_" + pv.dr_var + std::string(sign) + pair_suffix + "_platnorm" + cat;
                        gn->Write(save_name.c_str());

                        any_drawn = true;
                    }
                    if (!any_drawn) { delete c; continue; }

                    auto line = new TLine(pv.xmin, 1.0, pv.xmax, 1.0);
                    line->SetLineStyle(2);
                    line->SetLineColor(kGray+2);
                    line->Draw();
                    leg->Draw();

                    TLatex tex;
                    tex.SetNDC();
                    tex.SetTextSize(0.035);
                    std::string header = sign_label + ", pair-level, " + mode_label;
                    if (!cat_label.empty()) header += ", " + cat_label;
                    tex.DrawLatex(0.15, 0.88, header.c_str());

                    std::string outdir = base_outdir + "combined/" + ctr_subdir;
                    gSystem->mkdir(outdir.c_str(), true);
                    std::string outname = outdir + "/" + pv.dr_var + std::string(sign) + pair_suffix + "_platnorm.png";
                    c->SaveAs(outname.c_str());
                    delete c;
                }

                // --- Plot pair-pT-sliced DR_zoomin (normalized) ---
                for (auto& y : years) {
                    auto c = new TCanvas("c", "", 800, 600);
                    c->SetLeftMargin(0.12);
                    c->SetRightMargin(0.05);
                    auto frame = c->DrawFrame(0.0, 0.85, 0.8, 1.10);
                    frame->GetXaxis()->SetTitle("#DeltaR");
                    frame->GetYaxis()->SetTitle("#varepsilon_{#DeltaR}^{pair} (plateau normalized)");

                    auto leg = new TLegend(0.45, 0.15, 0.90, 0.45);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);

                    int colors[] = {kBlue, kRed, kGreen+2, kMagenta+1};
                    int markers[] = {20, 21, 22, 23};
                    bool any_drawn = false;

                    for (int isl = 0; isl < (int)pt_slices.size(); isl++) {
                        std::string gname = "g_DR_zoomin" + std::string(sign) + pair_suffix + "_invw_num" + cat + pt_slices[isl].suffix + "_divided";
                        auto g = (TGraphAsymmErrors*)files[y.year]->Get(gname.c_str());
                        if (!g) continue;

                        double plat;
                        if (use_pt_binned_plateau) {
                            plat = pt_bin_plateaus[y.year][pt_slices[isl].suffix];
                        } else {
                            plat = pt_int_plateau[y.year];
                        }
                        auto gn = normalizeGraph(g, plat);
                        gn->SetMarkerColor(colors[isl]);
                        gn->SetLineColor(colors[isl]);
                        gn->SetMarkerStyle(markers[isl]);
                        gn->SetMarkerSize(0.8);
                        gn->Draw("P same");
                        leg->AddEntry(gn, Form("%s (plat=%.3f)", pt_slices[isl].label.c_str(), plat), "lp");

                        out_files[y.year]->cd();
                        std::string save_name = "g_DR_zoomin" + std::string(sign) + pair_suffix + "_platnorm" + cat + pt_slices[isl].suffix;
                        gn->Write(save_name.c_str());

                        any_drawn = true;
                    }
                    if (!any_drawn) { delete c; continue; }

                    auto line = new TLine(0.0, 1.0, 0.8, 1.0);
                    line->SetLineStyle(2);
                    line->SetLineColor(kGray+2);
                    line->Draw();
                    leg->Draw();

                    TLatex tex;
                    tex.SetNDC();
                    tex.SetTextSize(0.035);
                    std::string header = Form("PbPb 20%d, %s, pair-level, %s", y.year, sign_label.c_str(), mode_label.c_str());
                    if (!cat_label.empty()) header += ", " + cat_label;
                    tex.DrawLatex(0.15, 0.88, header.c_str());

                    std::string outdir = base_outdir + "combined/" + (cat.empty() ? "ctr_integrated" : "ctr_binned");
                    gSystem->mkdir(outdir.c_str(), true);
                    std::string outname = outdir + "/DR_zoomin" + std::string(sign) + pair_suffix + "_platnorm" + cat + Form("_yr%d_pt_slices.png", y.year);
                    c->SaveAs(outname.c_str());
                    delete c;
                }
            }
        }

        for (auto& [yr, f] : out_files) {
            f->Write();
            f->Close();
            std::cout << "Written ROOT file: " << f->GetName() << std::endl;
        }
    };

    std::cout << "\n===== Mode 1: pT-integrated plateau normalization =====\n" << std::endl;
    makeNormPlots("pT-int plateau", "plateau_norm_pt_int", false);

    std::cout << "\n===== Mode 2: pT-binned plateau normalization =====\n" << std::endl;
    makeNormPlots("pT-binned plateau", "plateau_norm_pt_binned", true);

    for (auto& [yr, f] : files) f->Close();
    std::cout << "\nAll plateau-normalized plots saved." << std::endl;
}
