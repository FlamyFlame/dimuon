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

bool isBNL = true;

std::string bnl_base() { return "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/"; }
std::string mac_base() { return "/Users/yuhanguo/Documents/physics/heavy-ion/dimuon/datasets/"; }
std::string base_dir() { return isBNL ? bnl_base() : mac_base(); }

std::string plot_base() {
    return isBNL ? bnl_base() + "plots/" : mac_base() + "plots/";
}

void plot_dR_trig_corr() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    struct YearInfo {
        int year;
        std::string path;
        int color;
        int marker;
    };

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

    struct PlotDef {
        std::string var;
        std::string xtitle;
        double xmin, xmax;
        bool logx;
    };

    std::vector<PlotDef> plots = {
        {"DR",           "#DeltaR",              0.0, 5.75, false},
        {"DR_zoomin",    "#DeltaR",              0.0, 0.8,  false},
        {"DR_0_2",       "#DeltaR",              0.0, 0.2,  false},
        {"Deta_zoomin",  "#Delta#eta",          -0.8, 0.8,  false},
        {"Dphi_zoomin",  "#Delta#phi",          -0.8, 0.8,  false},
        {"pair_pt_log",  "p_{T}^{pair} [GeV]",  4.0, 80,   true},
        {"minv_zoomin",  "m_{#mu#mu} [GeV]",    0.3, 3.0,  false}
    };

    struct TermDef {
        std::string graph_suffix;
        std::string label;
        std::string subdir;
        double ymin, ymax;
    };

    std::vector<TermDef> terms = {
        {"_mu4_combined_invw_num",   "single-muon (combined)", "dR_single_muon", 0.6, 1.6},
        {"_cross_mu4_invw_num",      "cross-term",             "dR_cross_term",  0.5, 2.5}
    };

    auto makePlots = [&](const std::string& ctr_label, const std::string& ctr_graph_suffix,
                         const std::string& ctr_subdir) {
        for (const auto& sign : {"_op", "_ss"}) {
            std::string sign_label = (std::string(sign) == "_op") ? "OS" : "SS";

            for (const auto& pd : plots) {
                for (const auto& term : terms) {
                    std::string gname = "g_" + pd.var + sign + term.graph_suffix + ctr_graph_suffix + "_divided";

                    auto c = new TCanvas("c", "", 800, 600);
                    c->SetLeftMargin(0.12);
                    c->SetRightMargin(0.05);
                    if (pd.logx) c->SetLogx();

                    auto frame = c->DrawFrame(pd.xmin, term.ymin, pd.xmax, term.ymax);
                    frame->GetXaxis()->SetTitle(pd.xtitle.c_str());
                    frame->GetYaxis()->SetTitle("#varepsilon_{#DeltaR} (inv-weighted ratio)");

                    auto leg = new TLegend(0.55, 0.72, 0.90, 0.90);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);

                    bool any_drawn = false;
                    for (auto& y : years) {
                        auto g = (TGraphAsymmErrors*)files[y.year]->Get(gname.c_str());
                        if (!g) {
                            std::cerr << "WARNING: missing " << gname << " in year " << y.year << std::endl;
                            continue;
                        }
                        g->SetMarkerColor(y.color);
                        g->SetLineColor(y.color);
                        g->SetMarkerStyle(y.marker);
                        g->SetMarkerSize(0.8);
                        g->Draw("P same");
                        leg->AddEntry(g, Form("PbPb 20%d", y.year), "lp");
                        any_drawn = true;
                    }

                    if (!any_drawn) { delete c; continue; }

                    auto line = new TLine(pd.xmin, 1.0, pd.xmax, 1.0);
                    line->SetLineStyle(2);
                    line->SetLineColor(kGray+2);
                    line->Draw();

                    leg->Draw();

                    TLatex tex;
                    tex.SetNDC();
                    tex.SetTextSize(0.035);
                    std::string header = sign_label + ", " + term.label;
                    if (!ctr_label.empty()) header += ", " + ctr_label;
                    tex.DrawLatex(0.15, 0.88, header.c_str());

                    std::string outdir = plot_base() + "pbpb_trigger_efficiency/mu4/" + term.subdir + "/combined/" + ctr_subdir;
                    gSystem->mkdir(outdir.c_str(), true);
                    std::string outname = outdir + "/" + pd.var + sign + term.graph_suffix + ctr_graph_suffix + ".png";
                    c->SaveAs(outname.c_str());
                    delete c;
                }
            }
        }
    };

    makePlots("", "", "ctr_integrated");

    std::vector<std::pair<std::string, std::string>> ctr_bins = {
        {"0-5%",   "_ctr0_5"},
        {"5-10%",  "_ctr5_10"},
        {"10-20%", "_ctr10_20"},
        {"20-30%", "_ctr20_30"},
        {"30-50%", "_ctr30_50"},
        {"50-80%", "_ctr50_80"}
    };
    for (const auto& [label, suffix] : ctr_bins) {
        makePlots(label, suffix, "ctr_binned");
    }

    for (auto& [yr, f] : files) f->Close();
    std::cout << "All dR correction plots saved under " << plot_base() << "pbpb_trigger_efficiency/mu4/" << std::endl;
}
