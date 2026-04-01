// RecoEffyRetRespPlotter.h (or put in a .cxx you .L in ROOT)
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TError.h>
#include <TLatex.h>
#include <TPad.h>
#include <TGraphAsymmErrors.h>
#include <TRatioPlot.h>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <iostream>
#include "Utilities/PlotUtils.h"
#include "Utilities/PlotCommonConfig.h"

class RecoEffyRetRespPlotter {
public:

    RecoEffyRetRespPlotter(int run_year_input, bool tight_WP_input = false, bool require_signal_cuts_input = false)
    :   run_year(run_year_input % 2000),
        tight_WP(tight_WP_input),
        require_signal_cuts(require_signal_cuts_input)
    {
        isRun3 = !(run_year <= 18);
    }

    RecoEffyRetRespPlotter()
    : RecoEffyRetRespPlotter(17, false, false){}

    ~RecoEffyRetRespPlotter() {
        if (infile_unmixed) {
            infile_unmixed->Close();
            delete infile_unmixed;
            infile_unmixed = nullptr;
        }
        if (infile_mixed) {
            infile_mixed->Close();
            delete infile_mixed;
            infile_mixed = nullptr;
        }
    }

    void Run(){
        if (!Initialize()) return;

        PlotTruthRecoCompr();
        PlotResponseMatrix();
        Plot1DRecoEffcy();
        Plot2DRecoEffcySigned();
        Plot2DRecoEffcySingleB();
        Plot1DRecoEffcySingleBOpCompr();
        Plot1DRecoEffcyRangedSingleBOpCompr();
        if (!require_signal_cuts) Plot1DRecoEffcyRangedMixed();
    }

protected:
    int run_year = 17;
    bool isRun3;
    bool tight_WP;
    bool require_signal_cuts;

    std::string run_period;
    std::string wp_suffix; // working point suffix
    std::string require_signal_cuts_hist_suffix_num;
    std::string require_signal_cuts_hist_suffix_denom;
    std::string require_signal_cuts_file_dir_suffix;
    std::string wp_filter; // working point suffix

    PlotCommonConfig cfg;

    std::string data_dir;
    std::string infile_name_unmixed;
    std::string infile_name_mixed;
    TFile* infile_unmixed{nullptr};
    TFile* infile_mixed{nullptr};
    std::vector<std::string> pair_vars_for_det_resp = {"pair_pt", "minv_zoomin", "dr_zoomin"};
    // std::vector<std::string> pair_vars_for_1d_reco_effcy = {"truth_pair_pt", "truth_pair_eta", "truth_dr_zoomin"};
    std::vector<std::string> pair_vars_for_1d_reco_effcy = {"truth_pair_pt", "truth_pair_eta", "truth_dr_zoomin", "truth_dr_2_0", "truth_minv_zoomin"};

    std::vector<std::array<std::string,2>> pair_vars_for_2d_reco_effcy = {
        {"truth_pair_pt", "truth_pair_eta"}, 
        {"truth_pair_pt", "truth_dr_zoomin"}, 
        {"truth_pair_eta", "truth_dr_zoomin"},
        {"truth_deta_zoomin", "truth_dphi_zoomin"},
        {"truth_pair_pt", "truth_minv_zoomin"}, 
        {"truth_minv_zoomin", "truth_dr_zoomin"}
    };

    std::vector<std::pair<float, float>> dr_ranges_for_reco_effcy = {
        {0.0f, 0.2f},
        {0.2f, 0.4f},
        {0.4f, 0.6f},
        {0.6f, 1.0f}
    };

    std::vector<std::pair<float, float>> pair_pT_ranges_for_reco_effcy_dR = {
        {8.0f, 12.0f},
        {12.0f, 20.0f},
        {20.0f, std::numeric_limits<float>::max()}
    };

    std::vector<std::pair<float, float>> dr_ranges_for_reco_effcy_mixed_wide = {
        {0.0f, 1.0f},
        {1.0f, 2.5f},
        {2.5f, 4.0f}
    };

    std::vector<std::tuple<std::string, std::string, bool, const std::vector<std::pair<float, float>>*>>
        reco_eff_proj_divide_cfgs = {
            {"truth_pair_pt",  "truth_dr_zoomin", false, &dr_ranges_for_reco_effcy},
            {"truth_pair_eta", "truth_dr_zoomin", false, &dr_ranges_for_reco_effcy},
            {"truth_pair_pt",  "truth_dr_zoomin", true,  &pair_pT_ranges_for_reco_effcy_dR}
        };
    

    using QEtaBinning = std::vector<std::pair<float, float>>;
    QEtaBinning q_eta_proj_ranges_fine_excl_gap_run2 = {
        {-2.4f, -2.0f}, 
        {-2.0f, -1.6f}, 
        {-1.6f, -1.3f}, 
        {-0.9f, -0.5f}, 
        {-0.5f, -0.1f}, 
        {0.1f, 0.5f}, 
        {0.5f, 1.0f}, 
        {1.3f, 1.6f}, 
        {1.6f, 2.0f}, 
        {2.0f, 2.4f}
    };

    bool Initialize() {
        run_period = isRun3? "run3" : "run2";

        require_signal_cuts_hist_suffix_num = require_signal_cuts ? "_pass_signal_truth" : "";
        require_signal_cuts_hist_suffix_denom = require_signal_cuts ? "_and_signal_truth_and_reco" : "";
        require_signal_cuts_file_dir_suffix = require_signal_cuts ? "_require_signal_cuts" : "";

        wp_suffix = tight_WP? "_tightWP" : "_mediumWP";
        wp_filter = tight_WP? "_pass_tight" : "_pass_medium";

        data_dir = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample";
        infile_name_unmixed = "histograms_powheg_fullsim_pp" + std::to_string(run_year) + ".root";
        infile_name_mixed = "histograms_powheg_fullsim_pp" + std::to_string(run_year) + "_mixed.root";

        std::string infile_path_unmixed = data_dir;
        if (!infile_path_unmixed.empty() && infile_path_unmixed.back() != '/') infile_path_unmixed += '/';
        infile_path_unmixed += infile_name_unmixed;

        std::string infile_path_mixed = data_dir;
        if (!infile_path_mixed.empty() && infile_path_mixed.back() != '/') infile_path_mixed += '/';
        infile_path_mixed += infile_name_mixed;

        infile_unmixed = TFile::Open(infile_path_unmixed.c_str(), "READ");
        if (!infile_unmixed || infile_unmixed->IsZombie()) {
            std::cerr << "[RecoEffyRetRespPlotter::Initialize] ERROR: failed to open unmixed file: "
                      << infile_path_unmixed << "\n";
            if (infile_unmixed) { delete infile_unmixed; infile_unmixed = nullptr; }
            return false;
        }

        infile_mixed = TFile::Open(infile_path_mixed.c_str(), "READ");
        if (!infile_mixed || infile_mixed->IsZombie()) {
            std::cerr << "[RecoEffyRetRespPlotter::Initialize] ERROR: failed to open mixed file: "
                      << infile_path_mixed << "\n";
            if (infile_mixed) { delete infile_mixed; infile_mixed = nullptr; }
            return false;
        }

        return true;
    }

    void PlotTruthRecoCompr() {
        if (!CheckUnmixedFile()) return;

        const std::string outdir = GetDetRespOutDir();
        gSystem->mkdir(outdir.c_str(), kTRUE);

        for (const auto& var : pair_vars_for_det_resp) {
            const std::string h_truth_name = "h_truth_" + var + "_single_b" + wp_filter;
            const std::string h_reco_name  = "h_"       + var + "_single_b" + wp_filter;

            TH1* h_truth_in = GetHistUnmixed<TH1>(h_truth_name);
            TH1* h_reco_in  = GetHistUnmixed<TH1>(h_reco_name);

            if (!h_truth_in) {
                std::cerr << "[PlotTruthRecoCompr] WARNING: missing " << h_truth_name << "\n";
                continue;
            }
            if (!h_reco_in) {
                std::cerr << "[PlotTruthRecoCompr] WARNING: missing " << h_reco_name << "\n";
                continue;
            }

            // Work with detached clones so styling/division does not mutate input-file objects.
            TH1* h_truth = dynamic_cast<TH1*>(h_truth_in->Clone((h_truth_name + "_clone").c_str()));
            TH1* h_reco  = dynamic_cast<TH1*>(h_reco_in->Clone((h_reco_name + "_clone").c_str()));
            if (!h_truth || !h_reco) {
                std::cerr << "[PlotTruthRecoCompr] WARNING: failed to clone input histograms for " << var << "\n";
                delete h_truth;
                delete h_reco;
                continue;
            }
            h_truth->SetDirectory(nullptr);
            h_reco->SetDirectory(nullptr);
            if (h_truth->GetSumw2N() == 0) h_truth->Sumw2();
            if (h_reco->GetSumw2N() == 0) h_reco->Sumw2();

            TCanvas c("c", "c", 900, 760);
            c.cd();

            bool logy = LogAx(var, cfg);

            h_truth->SetLineWidth(2);
            h_truth->SetLineColor(kBlack);
            h_truth->SetStats(0);
            h_truth->SetTitle("");

            h_reco->SetLineWidth(2);
            h_reco->SetLineColor(kRed);
            h_reco->SetStats(0);
            h_reco->SetTitle("");

            auto min_positive_bin_content = [](const TH1* h) {
                double vmin = std::numeric_limits<double>::max();
                if (!h) return vmin;
                for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
                    const double v = h->GetBinContent(ib);
                    if (v > 0.0 && v < vmin) vmin = v;
                }
                return vmin;
            };

            // Build ratio as h_reco/h_truth using TH1::Divide semantics (symmetric errors).
            TRatioPlot rp(h_reco, h_truth, "divsym");
            rp.SetLeftMargin(0.12);
            rp.SetRightMargin(0.05);
            rp.SetUpTopMargin(0.08);
            rp.SetUpBottomMargin(0.02);
            rp.SetLowTopMargin(0.02);
            rp.SetLowBottomMargin(0.32);
            rp.SetH1DrawOpt("hist");
            rp.SetH2DrawOpt("hist");
            rp.SetGraphDrawOpt("PE1"); // points with error bars in ratio pad
            rp.Draw();
            rp.SetSplitFraction(0.30);
            rp.SetSeparationMargin(0.02);

            if (rp.GetUpperPad()) {
                rp.GetUpperPad()->SetTicks(1, 1);
                rp.GetUpperPad()->SetLogy(logy);
                rp.GetUpperPad()->cd();
            }

            const double maxy = std::max(h_truth->GetMaximum(), h_reco->GetMaximum());
            if (maxy > 0.0 && rp.GetUpperRefYaxis()) {
                rp.GetUpperRefYaxis()->SetTitleSize(0.050);
                rp.GetUpperRefYaxis()->SetTitleSize(0.038);
                rp.GetUpperRefYaxis()->SetLabelSize(0.030);
                rp.GetUpperRefYaxis()->SetTitleOffset(1.15);
            }
            if (rp.GetUpperRefXaxis()) {
                rp.GetUpperRefXaxis()->SetTitleSize(0.0);
                rp.GetUpperRefXaxis()->SetLabelSize(0.0);
            }

            if (rp.GetLowerRefYaxis()) {
                rp.GetLowerRefYaxis()->SetTitle("Reco/Truth");
                if (var == "minv_zoomin" || var == "dr_zoomin") rp.GetLowerRefYaxis()->SetNdivisions(405);
                else rp.GetLowerRefYaxis()->SetNdivisions(505);
                rp.GetLowerRefYaxis()->SetTitleSize(0.050);
                rp.GetLowerRefYaxis()->SetTitleOffset(0.95);
                rp.GetLowerRefYaxis()->SetLabelSize(0.040);
            }
            // Keep ROOT's automatic ratio range to avoid internal axis-bound warnings.
            if (rp.GetLowerRefXaxis()) {
                rp.GetLowerRefXaxis()->SetTitle(h_truth->GetXaxis()->GetTitle());
                rp.GetLowerRefXaxis()->SetTitleSize(0.060);
                rp.GetLowerRefXaxis()->SetTitleOffset(1.10);
                rp.GetLowerRefXaxis()->SetLabelSize(0.050);
            }

            if (rp.GetUpperPad()) {
                rp.GetUpperPad()->cd();
            }

            TLegend* leg = new TLegend(0.62, 0.7, 0.88, 0.85);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(h_truth, "Truth", "l");
            leg->AddEntry(h_reco,  "Reco",  "l");
            leg->Draw("same");

            TLatex lab;
            lab.SetNDC();
            lab.SetTextSize(0.045);
            lab.DrawLatex(0.12, 0.92, (var + " truth vs reco").c_str());

            c.Modified();
            c.Update();

            const std::string out = outdir + var + "_truth_reco_compr" + wp_suffix + ".png";
            c.SaveAs(out.c_str());

            delete h_truth;
            delete h_reco;
        }
    }

    void PlotResponseMatrix() {
        if (!CheckUnmixedFile()) return;

        const std::string outdir = GetDetRespOutDir();
        gSystem->mkdir(outdir.c_str(), kTRUE);

        for (const auto& var : pair_vars_for_det_resp) {
            const std::string h2_name =
                "h_" + var + "_vs_truth_" + var + "_single_b" + wp_filter;

            TH2* h2 = GetHistUnmixed<TH2>(h2_name);
            if (!h2) {
                std::cerr << "[PlotResponseMatrix] WARNING: missing " << h2_name << "\n";
                continue;
            }

            TCanvas c("c", "c", 850, 700);
            c.SetRightMargin(0.14);
            c.SetTicks(1,1);
            c.cd();
            gPad->SetLogz(true);
            if (var == "pair_pt") {
                gPad->SetLogx(true);
                gPad->SetLogy(true);
            }

            h2->SetTitle((var + " response matrix").c_str());
            h2->SetStats(0);
            h2->Draw("colz");

            const std::string out = outdir + var + "_response_matrix" + wp_suffix + ".png";
            c.SaveAs(out.c_str());
        }
    }

    void Plot1DRecoEffcy()
    {
        if (!CheckRecoFiles()) return;

        const std::string outdir = GetRecoEffcyOutDir();
        gSystem->mkdir(outdir.c_str(), kTRUE);

        for (const auto& var : pair_vars_for_1d_reco_effcy) {

            TCanvas c(Form("c_reco_eff1d_%s", var.c_str()), Form("c_reco_eff1d_%s", var.c_str()), 1400, 500);
            c.Divide(2,1);

            // ---- LEFT: same sign ----
            c.cd(1);
            gPad->SetTicks(1, 1);

            bool logx = LogAx(var, cfg);
            gPad->SetLogx(logx);

            const std::string num_filter = RecoEffNumFilter();
            const std::string denom_filter = RecoEffDenomFilter();

            const std::string hnum_name_ss   = "h_" + var + "_ss" + num_filter;
            const std::string hdenom_name_ss = "h_" + var + "_ss" + denom_filter;

            TH1D* h_num_ss   = GetHistReco<TH1D>(hnum_name_ss);
            TH1D* h_den_ss   = GetHistReco<TH1D>(hdenom_name_ss);

            if (!h_num_ss) {
                std::cerr << "[Plot1DRecoEffcy] WARNING: missing " << hnum_name_ss << "\n";
            }
            if (!h_den_ss) {
                std::cerr << "[Plot1DRecoEffcy] WARNING: missing " << hdenom_name_ss << "\n";
            }

            TH1D* h_eff_ss = nullptr;
            if (h_num_ss && h_den_ss) {
                h_eff_ss = dynamic_cast<TH1D*>(h_num_ss->Clone(("h_eff_ss_" + var).c_str()));
                h_eff_ss->SetDirectory(nullptr);
                h_eff_ss->SetStats(0);
                h_eff_ss->Divide(h_den_ss);

                h_eff_ss->SetTitle("same sign");
                h_eff_ss->GetXaxis()->SetTitle(h_den_ss->GetXaxis()->GetTitle());
                h_eff_ss->GetYaxis()->SetTitle("#varepsilon");
                h_eff_ss->SetMinimum(0.0);
                h_eff_ss->SetMaximum(1.0);

                h_eff_ss->SetLineColor(kBlack);
                h_eff_ss->SetMarkerColor(kBlack);
                h_eff_ss->SetLineWidth(2);
                h_eff_ss->SetMarkerStyle(20);
                h_eff_ss->Draw("E1");

                // // subplot title
                // TLatex lab;
                // lab.SetNDC();
                // lab.SetTextSize(0.045);
                // lab.DrawLatex(0.12, 0.92, "same sign");

                if (logx) AdjustLogXRangeForHist(h_eff_ss, h_den_ss);
            }

            // ---- RIGHT: opposite sign ----
            c.cd(2);
            gPad->SetTicks(1, 1);
            gPad->SetLogx(logx);

            const std::string hnum_name_op   = "h_" + var + "_op" + num_filter;
            const std::string hdenom_name_op = "h_" + var + "_op" + denom_filter;

            TH1D* h_num_op = GetHistReco<TH1D>(hnum_name_op);
            TH1D* h_den_op = GetHistReco<TH1D>(hdenom_name_op);

            if (!h_num_op) {
                std::cerr << "[Plot1DRecoEffcy] WARNING: missing " << hnum_name_op << "\n";
            }
            if (!h_den_op) {
                std::cerr << "[Plot1DRecoEffcy] WARNING: missing " << hdenom_name_op << "\n";
            }

            TH1D* h_eff_op = nullptr;
            if (h_num_op && h_den_op) {
                h_eff_op = dynamic_cast<TH1D*>(h_num_op->Clone(("h_eff_op_" + var).c_str()));
                h_eff_op->SetDirectory(nullptr);
                h_eff_op->SetStats(0);
                h_eff_op->Divide(h_den_op);

                h_eff_op->SetTitle("opposite sign");
                h_eff_op->GetXaxis()->SetTitle(h_den_op->GetXaxis()->GetTitle());
                h_eff_op->GetYaxis()->SetTitle("#varepsilon");
                h_eff_op->SetMinimum(0.0);
                h_eff_op->SetMaximum(1.0);

                h_eff_op->SetLineColor(kBlack);
                h_eff_op->SetMarkerColor(kBlack);
                h_eff_op->SetLineWidth(2);
                h_eff_op->SetMarkerStyle(20);
                h_eff_op->Draw("E1");

                TLatex lab;
                lab.SetNDC();
                lab.SetTextSize(0.045);
                lab.DrawLatex(0.12, 0.92, "opposite sign");

                if (logx) AdjustLogXRangeForHist(h_eff_op, h_den_op);
            }

            // save
            const std::string out = outdir + "reco_effcy_" + var + wp_suffix + ".png";
            c.SaveAs(out.c_str());

            delete h_eff_ss;
            delete h_eff_op;
        }
    }

    void Plot2DRecoEffcySigned()
    {
        if (!CheckRecoFiles()) return;

        const std::string outdir = GetRecoEffcyOutDir();
        gSystem->mkdir(outdir.c_str(), kTRUE);

        for (const auto& var_pair : pair_vars_for_2d_reco_effcy) {
            const std::string& varx = var_pair[0];
            const std::string& vary = var_pair[1];

            TCanvas c(Form("c_reco_eff2d_%s_vs_%s", vary.c_str(), varx.c_str()), Form("c_reco_eff2d_%s_vs_%s", vary.c_str(), varx.c_str()), 1300, 550);
            c.Divide(2, 1);

            // "similar to Plot1DRecoEffcy": use logx decision via map; also allow logy for y-var if present
            const bool logx = LogAx(varx, cfg);
            const bool logy = LogAx(vary, cfg);
            const std::string num_filter = RecoEffNumFilter();
            const std::string denom_filter = RecoEffDenomFilter();

            // ---------- left pad: same sign ----------
            c.cd(1);
            gPad->SetTicks(1, 1);
            gPad->SetLogx(logx);
            gPad->SetLogy(logy);
            gPad->SetRightMargin(0.14);

            const std::string hnum_ss   = "h_" + vary + "_vs_" + varx + "_ss" + num_filter;
            const std::string hden_ss   = "h_" + vary + "_vs_" + varx + "_ss" + denom_filter;

            TH2D* h_num_ss = GetHistReco<TH2D>(hnum_ss);
            TH2D* h_den_ss = GetHistReco<TH2D>(hden_ss);

            if (!h_num_ss) { std::cerr << "[Plot2DRecoEffcy] WARNING: missing " << hnum_ss << "\n"; }
            if (!h_den_ss) { std::cerr << "[Plot2DRecoEffcy] WARNING: missing " << hden_ss << "\n"; }

            TH2D* h_eff_ss = nullptr;
            if (h_num_ss && h_den_ss) {
                h_eff_ss = dynamic_cast<TH2D*>(h_num_ss->Clone(("h_eff_ss_" + vary + "_vs_" + varx).c_str()));
                h_eff_ss->SetDirectory(nullptr); // avoid being tied to infile
                h_eff_ss->SetStats(0);
                h_eff_ss->Divide(h_den_ss);

                h_eff_ss->SetMinimum(0.0);
                h_eff_ss->SetMaximum(1.0);

                h_eff_ss->SetTitle("same sign");

                h_eff_ss->Draw("colz");

                // TLatex lab;
                // lab.SetNDC();
                // lab.SetTextSize(0.045);
                // lab.DrawLatex(0.12, 0.92, "same sign");
            }

            // ---------- right pad: opposite sign ----------
            c.cd(2);
            gPad->SetTicks(1, 1);
            gPad->SetLogx(logx);
            gPad->SetLogy(logy);
            gPad->SetRightMargin(0.14);

            const std::string hnum_op   = "h_" + vary + "_vs_" + varx + "_op" + num_filter;
            const std::string hden_op   = "h_" + vary + "_vs_" + varx + "_op" + denom_filter;

            TH2D* h_num_op = GetHistReco<TH2D>(hnum_op);
            TH2D* h_den_op = GetHistReco<TH2D>(hden_op);

            if (!h_num_op) { std::cerr << "[Plot2DRecoEffcy] WARNING: missing " << hnum_op << "\n"; }
            if (!h_den_op) { std::cerr << "[Plot2DRecoEffcy] WARNING: missing " << hden_op << "\n"; }

            TH2D* h_eff_op = nullptr;
            if (h_num_op && h_den_op) {
                h_eff_op = dynamic_cast<TH2D*>(h_num_op->Clone(("h_eff_op_" + vary + "_vs_" + varx).c_str()));
                h_eff_op->SetDirectory(nullptr);
                h_eff_op->SetStats(0);
                h_eff_op->Divide(h_den_op);

                h_eff_op->SetMinimum(0.0);
                h_eff_op->SetMaximum(1.0);

                h_eff_op->SetTitle("opposite sign");

                h_eff_op->Draw("colz");

                // TLatex lab;
                // lab.SetNDC();
                // lab.SetTextSize(0.045);
                // lab.DrawLatex(0.12, 0.92, "opposite sign");
            }

            const std::string out = outdir + "reco_effcy_" + vary + "_vs_" + varx + wp_suffix + ".png";
            c.SaveAs(out.c_str());

            delete h_eff_ss;
            delete h_eff_op;
        }
    }

    void Plot2DRecoEffcySingleB()
    {
        if (!CheckRecoFiles()) return;

        const std::string outdir = GetRecoEffcyOutDir();
        gSystem->mkdir(outdir.c_str(), kTRUE);

        for (const auto& var_pair : pair_vars_for_2d_reco_effcy) {
            const std::string& varx = var_pair[0];
            const std::string& vary = var_pair[1];

            TCanvas c(Form("c_reco_eff2d_%s_vs_%s_single_b", vary.c_str(), varx.c_str()),
                      Form("c_reco_eff2d_%s_vs_%s_single_b", vary.c_str(), varx.c_str()),
                      700, 600);
            c.cd();

            const bool logx = LogAx(varx, cfg);
            const bool logy = LogAx(vary, cfg);
            const std::string num_filter = RecoEffNumFilter();
            const std::string denom_filter = RecoEffDenomFilter();

            gPad->SetTicks(1, 1);
            gPad->SetLogx(logx);
            gPad->SetLogy(logy);
            gPad->SetRightMargin(0.14);

            const std::string hnum_single_b = "h_" + vary + "_vs_" + varx + "_single_b" + num_filter;
            const std::string hden_single_b = "h_" + vary + "_vs_" + varx + "_single_b" + denom_filter;

            TH2D* h_num_single_b = GetHistReco<TH2D>(hnum_single_b);
            TH2D* h_den_single_b = GetHistReco<TH2D>(hden_single_b);

            if (!h_num_single_b) { std::cerr << "[Plot2DRecoEffcySingleB] WARNING: missing " << hnum_single_b << "\n"; }
            if (!h_den_single_b) { std::cerr << "[Plot2DRecoEffcySingleB] WARNING: missing " << hden_single_b << "\n"; }

            TH2D* h_eff_single_b = nullptr;
            if (h_num_single_b && h_den_single_b) {
                h_eff_single_b = dynamic_cast<TH2D*>(h_num_single_b->Clone(("h_eff_single_b_" + vary + "_vs_" + varx).c_str()));
                h_eff_single_b->SetDirectory(nullptr);
                h_eff_single_b->SetStats(0);
                h_eff_single_b->Divide(h_den_single_b);

                h_eff_single_b->SetMinimum(0.0);
                h_eff_single_b->SetMaximum(1.0);
                h_eff_single_b->SetTitle("single-b");

                h_eff_single_b->Draw("colz");
            }

            const std::string out = outdir + "reco_effcy_" + vary + "_vs_" + varx + "_single_b" + wp_suffix + ".png";
            c.SaveAs(out.c_str());

            delete h_eff_single_b;
        }
    }

// Compare opposite-sign vs single-b signal reco efficiency (truth-level denominators)
// Style intentionally mirrors Plot1DRecoEffcy (ticks/logx/axis ranges).
void Plot1DRecoEffcySingleBOpCompr()
{
    if (!CheckRecoFiles()) return;

    const std::string outdir = GetRecoEffcyOutDir();
    gSystem->mkdir(outdir.c_str(), kTRUE);

    for (const auto& var : pair_vars_for_1d_reco_effcy) {

        TCanvas c(Form("c_reco_eff1d_%s_single_b_op_compr", var.c_str()),
                  Form("c_reco_eff1d_%s_single_b_op_compr", var.c_str()),
                  800, 600);
        c.cd();
        gPad->SetTicks(1, 1);

        const bool logx = LogAx(var, cfg);
        gPad->SetLogx(logx);

        const std::string num_filter = RecoEffNumFilter();
        const std::string denom_filter = RecoEffDenomFilter();

        // ---- opposite sign ----
        const std::string hnum_name_op   = "h_" + var + "_op" + num_filter;
        const std::string hdenom_name_op = "h_" + var + "_op" + denom_filter;

        TH1D* h_num_op = GetHistReco<TH1D>(hnum_name_op);
        TH1D* h_den_op = GetHistReco<TH1D>(hdenom_name_op);

        if (!h_num_op)  { std::cerr << "[Plot1DRecoEffcySingleBOpCompr] WARNING: missing " << hnum_name_op   << "\n"; }
        if (!h_den_op)  { std::cerr << "[Plot1DRecoEffcySingleBOpCompr] WARNING: missing " << hdenom_name_op << "\n"; }

        // ---- single-b ----
        const std::string hnum_name_single_b   = "h_" + var + "_single_b" + num_filter;
        const std::string hdenom_name_single_b = "h_" + var + "_single_b" + denom_filter;

        TH1D* h_num_single_b = GetHistReco<TH1D>(hnum_name_single_b);
        TH1D* h_den_single_b = GetHistReco<TH1D>(hdenom_name_single_b);

        if (!h_num_single_b) { std::cerr << "[Plot1DRecoEffcySingleBOpCompr] WARNING: missing " << hnum_name_single_b   << "\n"; }
        if (!h_den_single_b) { std::cerr << "[Plot1DRecoEffcySingleBOpCompr] WARNING: missing " << hdenom_name_single_b << "\n"; }

        TH1D* h_eff_op = nullptr;
        TH1D* h_eff_single_b = nullptr;

        if (h_num_op && h_den_op) {
            h_eff_op = dynamic_cast<TH1D*>(h_num_op->Clone(("h_eff_op_" + var).c_str()));
            h_eff_op->SetDirectory(nullptr);
            h_eff_op->SetStats(0);
            h_eff_op->Divide(h_den_op);

            h_eff_op->GetXaxis()->SetTitle(h_den_op->GetXaxis()->GetTitle());
            h_eff_op->SetLineColor(kBlack);
            h_eff_op->SetMarkerColor(kBlack);
            h_eff_op->SetLineWidth(2);
            h_eff_op->SetMarkerStyle(20);
        }
        if (h_num_single_b && h_den_single_b) {
            h_eff_single_b = dynamic_cast<TH1D*>(h_num_single_b->Clone(("h_eff_single_b_" + var).c_str()));
            h_eff_single_b->SetDirectory(nullptr);
            h_eff_single_b->SetStats(0);
            h_eff_single_b->Divide(h_den_single_b);

            h_eff_single_b->GetXaxis()->SetTitle(h_den_single_b->GetXaxis()->GetTitle());
            h_eff_single_b->SetLineColor(kRed);
            h_eff_single_b->SetMarkerColor(kRed);
            h_eff_single_b->SetLineWidth(2);
            h_eff_single_b->SetMarkerStyle(20);
        }

        if (h_eff_op || h_eff_single_b) {
            TH1D* h_frame = h_eff_op ? h_eff_op : h_eff_single_b;
            h_frame->GetYaxis()->SetTitle("#varepsilon");
            h_frame->SetMinimum(0.0);
            h_frame->SetMaximum(1.0);
            h_frame->Draw("E1");

            if (h_eff_op && h_eff_op != h_frame)                   h_eff_op->Draw("E1 same");
            if (h_eff_single_b && h_eff_single_b != h_frame)       h_eff_single_b->Draw("E1 same");

            // Match log-x range behavior used in Plot1DRecoEffcy.
            // Use opposite-sign denom as the reference; if missing, fall back to single-b denom.
            if (logx) {
                if (h_den_op) AdjustLogXRangeForHist(h_frame, h_den_op);
                else if (h_den_single_b) AdjustLogXRangeForHist(h_frame, h_den_single_b);
            }

            TLegend* leg = new TLegend(0.62, 0.7, 0.85, 0.85);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            if (h_eff_op)       leg->AddEntry(h_eff_op,       "opposite sign", "lp");
            if (h_eff_single_b) leg->AddEntry(h_eff_single_b, "single-b",      "lp");
            leg->Draw("same");
        }

        const std::string out = outdir + "reco_effcy_" + var + "_single_b_op_compr" + wp_suffix + ".png";
        c.SaveAs(out.c_str());

        delete h_eff_op;
        delete h_eff_single_b;
    }
}

void Plot1DRecoEffcyRangedSingleBOpComprHelper(
    const std::string& varx,
    const std::string& vary,
    bool proj_axis,
    const std::vector<std::pair<float, float>>& proj_ranges)
{
    if (!CheckRecoFiles()) return;
        const std::string num_filter = RecoEffNumFilter();

    if (proj_ranges.empty()) return;

    auto axisTag = [proj_axis]() {
        return proj_axis ? std::string("_py") : std::string("_px");
    };

    auto as_suffix = [](const std::pair<float, float>& range) {
        auto format_num = [](float x) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << x;
            std::string s = oss.str();
            for (auto& c : s) if (c == '.') c = '_';
            while (s.find('-') != std::string::npos) s.replace(s.find('-'), 1, "minus");
            return s;
        };

        const bool upper_is_max = (!std::isfinite(range.second)
                                || range.second >= std::numeric_limits<float>::max() * 0.5f);
        if (upper_is_max) return format_num(range.first) + "_TO_MAX";
        return format_num(range.first) + "_TO_" + format_num(range.second);
    };

    auto fmt_range_val = [](float x) {
        if (!std::isfinite(x) || x >= std::numeric_limits<float>::max() * 0.5f) return std::string("MAX");
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1) << x;
        std::string s = oss.str();
        if (s.find('.') != std::string::npos) {
            while (!s.empty() && s.back() == '0') s.pop_back();
            if (!s.empty() && s.back() == '.') s.pop_back();
        }
        return s;
    };

    auto var_label = [](const std::string& var) {
        if (var == "truth_dr_zoomin") return std::string("#Delta R");
        if (var == "truth_pair_pt") return std::string("p_{T}^{truth}");
        if (var == "truth_pair_eta") return std::string("#eta^{truth}");
        if (var == "truth_minv_zoomin") return std::string("m_{#mu#mu}^{truth}");
        return var;
    };

    auto best_grid = [](int n) {
        if (n <= 3) return std::pair<int,int>(1, n);

        int best_rows = n;
        int best_cols = 1;
        int best_diff = n - 1;
        int best_area = n;

        for (int cols = 1; cols <= n; ++cols) {
            int rows = (n + cols - 1) / cols;
            if (rows < cols) continue;

            int diff = rows - cols;
            int area = rows * cols;
            if (diff < best_diff || (diff == best_diff && area < best_area)) {
                best_rows = rows;
                best_cols = cols;
                best_diff = diff;
                best_area = area;
            }
        }
        return std::pair<int,int>(best_rows, best_cols);
    };

    const auto outdir = GetRecoEffcyOutDir();
    gSystem->mkdir(outdir.c_str(), kTRUE);
    const auto ranged_outdir = outdir + "ranged/";
    gSystem->mkdir(ranged_outdir.c_str(), kTRUE);

    const auto [nrows, ncols] = best_grid(static_cast<int>(proj_ranges.size()));
    TCanvas c(Form("c_reco_eff1d_ranged_%s_vs_%s%s", vary.c_str(), varx.c_str(), axisTag().c_str()),
              Form("c_reco_eff1d_ranged_%s_vs_%s%s", vary.c_str(), varx.c_str(), axisTag().c_str()),
              500 * ncols, 380 * nrows);
    c.Divide(ncols, nrows);

    const std::string xvar = proj_axis ? vary : varx;
    const bool projx = !proj_axis;
    const std::string var_target = projx ? varx : vary;
    const bool logx = LogAx(xvar, cfg);

    for (std::size_t i = 0; i < proj_ranges.size(); ++i) {
        c.cd(static_cast<int>(i) + 1);
        gPad->SetTicks(1, 1);
        gPad->SetLogx(logx);

        const auto& range = proj_ranges[i];
        const std::string range_suffix = as_suffix(range);

        const std::string hname_op = "h_" + vary + "_vs_" + varx + "_op" + num_filter + axisTag() + "_" + range_suffix + "_divided";
        const std::string hname_single_b = "h_" + vary + "_vs_" + varx + "_single_b" + num_filter + axisTag() + "_" + range_suffix + "_divided";

        TH1D* h_op = GetHistReco<TH1D>(hname_op);
        TH1D* h_single_b = GetHistReco<TH1D>(hname_single_b);

        if (!h_op) {
            std::cerr << "[Plot1DRecoEffcyRangedSingleBOpComprHelper] WARNING: missing " << hname_op << "\n";
        }
        if (!h_single_b) {
            std::cerr << "[Plot1DRecoEffcyRangedSingleBOpComprHelper] WARNING: missing " << hname_single_b << "\n";
        }

        if (!h_op && !h_single_b) continue;

        if (h_op) {
            h_op->SetLineColor(kBlack);
            h_op->SetMarkerColor(kBlack);
            h_op->SetLineWidth(2);
            h_op->SetMarkerStyle(20);
            h_op->SetTitle("");
            h_op->SetStats(0);
        }
        if (h_single_b) {
            h_single_b->SetLineColor(kRed);
            h_single_b->SetMarkerColor(kRed);
            h_single_b->SetLineWidth(2);
            h_single_b->SetMarkerStyle(20);
            h_single_b->SetTitle("");
            h_single_b->SetStats(0);
        }

        if (h_op) {
            h_op->GetYaxis()->SetTitle("#varepsilon");
            h_op->SetMinimum(0.0);
            h_op->SetMaximum(1.0);
            std::string x_title = h_op->GetXaxis()->GetTitle();
            if (x_title.empty()) x_title = var_target;
            h_op->GetXaxis()->SetTitle(x_title.c_str());
            h_op->Draw("E1");

            if (logx) {
                AdjustLogXRangeForHist(h_op, h_op);
            }
        } else if (h_single_b) {
            h_single_b->GetYaxis()->SetTitle("#varepsilon");
            h_single_b->SetMinimum(0.0);
            h_single_b->SetMaximum(1.0);
            std::string x_title = h_single_b->GetXaxis()->GetTitle();
            if (x_title.empty()) x_title = var_target;
            h_single_b->GetXaxis()->SetTitle(x_title.c_str());
            h_single_b->Draw("E1");

            if (logx) {
                AdjustLogXRangeForHist(h_single_b, h_single_b);
            }
        }

        if (h_single_b) h_single_b->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.12, 0.12, 0.54, 0.36);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.034);
        if (h_op) leg->AddEntry(h_op, "opposite sign", "lp");
        if (h_single_b) leg->AddEntry(h_single_b, "single-b", "lp");

        const std::string ranged_var = proj_axis ? varx : vary;
        const std::string ranged_unit = (ranged_var == "truth_pair_pt") ? " [GeV]" : "";
        const std::string ranged_text = var_label(ranged_var) + ": [" + fmt_range_val(range.first) + ", " + fmt_range_val(range.second) + "]" + ranged_unit;
        leg->AddEntry((TObject*)nullptr, ranged_text.c_str(), "");
        leg->Draw("same");
    }

    const std::string var_range = projx ? vary : varx;
    const std::string out = ranged_outdir + "reco_effcy_single_b_op_compr_" + var_target + "_in_" + var_range + "_bins" + wp_suffix + ".png";
    c.SaveAs(out.c_str());
}

void Plot1DRecoEffcyRangedSingleBOpCompr()
{
    if (!CheckRecoFiles()) return;

    for (const auto& cfg_tuple : reco_eff_proj_divide_cfgs) {
        const std::string& varx = std::get<0>(cfg_tuple);
        const std::string& vary = std::get<1>(cfg_tuple);
        const bool proj_axis = std::get<2>(cfg_tuple);
        const auto* proj_ranges = std::get<3>(cfg_tuple);

        if (proj_ranges == nullptr || proj_ranges->empty()) continue;
        Plot1DRecoEffcyRangedSingleBOpComprHelper(varx, vary, proj_axis, *proj_ranges);
    }
}

void Plot1DRecoEffcyRangedMixed()
{
    if (!CheckRecoFiles()) return;
    if (require_signal_cuts) return;

    const std::string outdir = GetRecoEffcyOutDir();
    gSystem->mkdir(outdir.c_str(), kTRUE);
    const std::string ranged_outdir = outdir + "ranged/";
    gSystem->mkdir(ranged_outdir.c_str(), kTRUE);

    auto as_suffix = [](const std::pair<float, float>& range) {
        auto format_num = [](float x) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << x;
            std::string s = oss.str();
            for (auto& c : s) if (c == '.') c = '_';
            while (s.find('-') != std::string::npos) s.replace(s.find('-'), 1, "minus");
            return s;
        };
        return format_num(range.first) + "_TO_" + format_num(range.second);
    };

    auto fmt_range = [](const std::pair<float, float>& range) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(1) << range.first << " < #DeltaR < " << range.second;
        return oss.str();
    };

    const std::vector<std::string> vars = {"truth_pair_pt", "truth_pair_eta"};
    const std::vector<std::string> signs = {"_ss", "_op"};
    const std::vector<Color_t> colors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 1};

    const std::string num_filter = RecoEffNumFilter(); // _pass_medium or _pass_tight

    for (const auto& varx : vars) {
        TCanvas c(Form("c_reco_eff_mixed_ranged_%s", varx.c_str()), Form("c_reco_eff_mixed_ranged_%s", varx.c_str()), 1400, 550);
        c.Divide(2, 1);

        for (std::size_t is = 0; is < signs.size(); ++is) {
            c.cd(static_cast<int>(is) + 1);
            gPad->SetTicks(1, 1);
            gPad->SetLogx(LogAx(varx, cfg));

            std::vector<TH1D*> hs;
            hs.reserve(dr_ranges_for_reco_effcy_mixed_wide.size());

            for (std::size_t ir = 0; ir < dr_ranges_for_reco_effcy_mixed_wide.size(); ++ir) {
                const auto& range = dr_ranges_for_reco_effcy_mixed_wide[ir];
                const std::string range_suffix = as_suffix(range);
                const std::string hname = "h_truth_dr_4_0_vs_" + varx + signs[is] + num_filter + "_px_" + range_suffix + "_divided";

                TH1D* h = GetHistReco<TH1D>(hname);
                if (!h) {
                    throw std::runtime_error("Plot1DRecoEffcyRangedMixed: missing histogram '" + hname + "' in mixed file.");
                }

                h->SetLineColor(colors[ir % colors.size()]);
                h->SetMarkerColor(colors[ir % colors.size()]);
                h->SetLineWidth(2);
                h->SetMarkerStyle(20);
                h->SetMarkerSize(0.8);
                h->SetStats(0);
                h->GetYaxis()->SetTitle("#varepsilon");
                h->SetMinimum(0.0);
                h->SetMaximum(1.0);
                hs.push_back(h);
            }

            if (!hs.empty()) {
                hs.front()->SetTitle((std::string(is == 0 ? "same-sign" : "opposite-sign") + ";" + hs.front()->GetXaxis()->GetTitle() + ";#varepsilon").c_str());
                hs.front()->Draw("E1");
                if (LogAx(varx, cfg)) AdjustLogXRangeForHist(hs.front(), hs.front());
                for (std::size_t ih = 1; ih < hs.size(); ++ih) hs[ih]->Draw("E1 SAME");

                TLegend* leg = new TLegend(0.14, 0.12, 0.58, 0.42);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.03);
                for (std::size_t ir = 0; ir < dr_ranges_for_reco_effcy_mixed_wide.size(); ++ir) {
                    leg->AddEntry(hs[ir], fmt_range(dr_ranges_for_reco_effcy_mixed_wide[ir]).c_str(), "lp");
                }
                leg->Draw("same");
            }
        }

        const std::string out = ranged_outdir + "reco_effcy_mixed_" + varx + "_in_truth_dr_bins" + wp_suffix + ".png";
        c.SaveAs(out.c_str());
    }
}

private:
    std::string RecoEffNumFilter() const {
        if (!require_signal_cuts) return wp_filter;
        return wp_filter + require_signal_cuts_hist_suffix_denom;
    }

    std::string RecoEffDenomFilter() const {
        if (!require_signal_cuts) return "";
        return require_signal_cuts_hist_suffix_num;
    }

    bool HistNameUsesMixed(const std::string& hname) const {
        const bool has_ss = (hname.find("_ss") != std::string::npos);
        const bool has_op = (hname.find("_op") != std::string::npos);
        const bool has_single_b = (hname.find("_single_b") != std::string::npos);

        if ((has_ss || has_op) && has_single_b) {
            throw std::runtime_error("[RecoEffyRetRespPlotter] Invalid reco histogram name with mixed+single_b tags: " + hname);
        }
        if (has_ss || has_op) return true;
        if (has_single_b) return false;

        throw std::runtime_error("[RecoEffyRetRespPlotter] Cannot determine reco input source for histogram name: " + hname);
    }

    template <typename T>
    T* GetHistUnmixed(const std::string& hname) const {
        if (!infile_unmixed) return nullptr;
        return dynamic_cast<T*>(infile_unmixed->Get(hname.c_str()));
    }

    template <typename T>
    T* GetHistReco(const std::string& hname) const {
        const bool use_mixed = HistNameUsesMixed(hname);
        TFile* src = use_mixed ? infile_mixed : infile_unmixed;
        if (!src) return nullptr;
        return dynamic_cast<T*>(src->Get(hname.c_str()));
    }

    void AdjustLogXRangeForHist(TH1* h_draw, const TH1* h_ref) const {
        if (!h_draw || !h_ref) return;

        const int nbins = h_ref->GetNbinsX();
        if (nbins <= 0) return;

        int first_bin = 1;
        while (first_bin <= nbins && h_ref->GetBinLowEdge(first_bin) <= 0.0) ++first_bin;
        if (first_bin > nbins) return;

        const double xmin = h_ref->GetBinLowEdge(first_bin);
        const double xmax = h_ref->GetBinLowEdge(nbins) + h_ref->GetBinWidth(nbins);
        h_draw->GetXaxis()->SetLimits(xmin, xmax);
    }

    std::string GetDetRespOutDir() const {
        // Ensure trailing slash on data_dir (optional)
        std::string base = data_dir;
        if (!base.empty() && base.back() != '/') base += '/';
        return base + "plots/" + run_period + "_det_resp_plots/" + (tight_WP ? "tight/" : "medium/");
    }

    bool CheckUnmixedFile() const {
        if (!infile_unmixed || !infile_unmixed->IsOpen()) {
            std::cerr << "[RecoEffyRetRespPlotter] ERROR: unmixed input file not open. Call Initialize().\n";
            return false;
        }
        return true;
    }

    bool CheckRecoFiles() const {
        if (!infile_unmixed || !infile_unmixed->IsOpen()) {
            std::cerr << "[RecoEffyRetRespPlotter] ERROR: unmixed input file not open. Call Initialize().\n";
            return false;
        }
        if (!infile_mixed || !infile_mixed->IsOpen()) {
            std::cerr << "[RecoEffyRetRespPlotter] ERROR: mixed input file not open. Call Initialize().\n";
            return false;
        }
        return true;
    }

    std::string GetRecoEffcyOutDir() const
    {
        std::string base = data_dir;
        if (!base.empty() && base.back() != '/') base += '/';
        return base + "plots/" + run_period + "_reco_effcy_plots" + require_signal_cuts_file_dir_suffix + "/" + (tight_WP ? "tight/" : "medium/");
    }

};
