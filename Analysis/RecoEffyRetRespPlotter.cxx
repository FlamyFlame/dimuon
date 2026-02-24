// DetRespPlotter.h (or put in a .cxx you .L in ROOT)
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TError.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <vector>
#include <iostream>
#include "Utilities/PlotUtils.h"
#include "Utilities/PlotCommonConfig.h"

class DetRespPlotter {
public:

    DetRespPlotter(int run_year_input, bool tight_WP_input = false, bool no_reco_resn_cuts_input = false)
    :   run_year(run_year_input % 2000),
        tight_WP(tight_WP_input),
        no_reco_resn_cuts(no_reco_resn_cuts_input)
    {
        isRun3 = !(run_year <= 18);
    }

    DetRespPlotter()
    : DetRespPlotter(17, false, false){}

    ~DetRespPlotter() {
        if (infile) {
            infile->Close();
            delete infile;
            infile = nullptr;
        }
    }

    void Run(){
        Initialize();
        if (!no_reco_resn_cuts){ // do not run detector response for pre-resonance cross-check          
            PlotTruthRecoCompr();
            PlotResponseMatrix();
        }
        Plot1DRecoEffcy();
        Plot2DRecoEffcySigned();
        Plot2DRecoEffcySingleB();
        Plot1DRecoEffcySingleBOpCompr();
        if (!no_reco_resn_cuts){
            Plot1DRecoEffcyRangedSingleBOpCompr();
        }
    }

protected:
    int run_year = 17;
    bool isRun3;
    bool tight_WP;
    bool no_reco_resn_cuts; // true if not apply reco resn cuts (cross-check for resn cuts' effects on reco effcies)

    std::string run_period;
    std::string wp_suffix; // working point suffix
    std::string no_reco_resn_cuts_hist_suffix; // histogram suffix if apply no reco resn cuts
    std::string no_reco_resn_cuts_file_dir_suffix; // file & dir suffix if apply no reco resn cuts
    std::string truth_resn_filter; // working point suffix
    std::string wp_filter; // working point suffix

    PlotCommonConfig cfg;

    std::string data_dir;
    std::string infile_name;
    TFile* infile{nullptr};
    std::vector<std::string> pair_vars_for_det_resp = {"pair_pt", "minv_zoomin", "dr_zoomin"};
    // std::vector<std::string> pair_vars_for_1d_reco_effcy = {"truth_pair_pt", "truth_pair_eta", "truth_dr_zoomin"};
    std::vector<std::string> pair_vars_for_1d_reco_effcy = {"truth_pair_pt", "truth_pair_eta", "truth_dr_zoomin", "truth_dr_2_0", "truth_minv_zoomin"};

    std::vector<std::array<std::string,2>> pair_vars_for_2d_reco_effcy = {
        {"truth_pair_pt", "truth_pair_eta"}, 
        {"truth_pair_pt", "truth_dr_zoomin"}, 
        {"truth_pair_eta", "truth_dr_zoomin"},
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

        no_reco_resn_cuts_hist_suffix = no_reco_resn_cuts ? "_pre_resonance" : "";
        no_reco_resn_cuts_file_dir_suffix = no_reco_resn_cuts ? "_no_reco_resn_cuts" : "";

        wp_suffix = tight_WP? "_tightWP" : "_mediumWP";
        wp_filter = tight_WP? "_pair_pass_tight" + no_reco_resn_cuts_hist_suffix 
                            : "_pair_pass_medium" + no_reco_resn_cuts_hist_suffix;
        truth_resn_filter = "_pair_pass_resonance_truth"; // set to "" to turn off

        data_dir = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample";
        infile_name = "histograms_powheg_fullsim_pp" + std::to_string(run_year) + no_reco_resn_cuts_file_dir_suffix + ".root";

        std::string infile_path = data_dir;
        if (!infile_path.empty() && infile_path.back() != '/') infile_path += '/';
        infile_path += infile_name;
        // Open input file
        infile = TFile::Open(infile_path.c_str(), "READ");
        if (!infile || infile->IsZombie()) {
            std::cerr << "[DetRespPlotter::Initialize] ERROR: failed to open file: "
                      << infile_path << "\n";
            if (infile) { delete infile; infile = nullptr; }
            return false;
        }

        return true;
    }

    void PlotTruthRecoCompr() {
        if (!CheckFile()) return;

        const std::string outdir = GetDetRespOutDir();
        makeDirIfNeeded(outdir);

        for (const auto& var : pair_vars_for_det_resp) {
            const std::string h_truth_name = "h_truth_" + var + "_single_b" + wp_filter;
            const std::string h_reco_name  = "h_"       + var + "_single_b" + wp_filter;

            TH1* h_truth = dynamic_cast<TH1*>(infile->Get(h_truth_name.c_str()));
            TH1* h_reco  = dynamic_cast<TH1*>(infile->Get(h_reco_name.c_str()));

            if (!h_truth) {
                std::cerr << "[PlotTruthRecoCompr] WARNING: missing " << h_truth_name << "\n";
                continue;
            }
            if (!h_reco) {
                std::cerr << "[PlotTruthRecoCompr] WARNING: missing " << h_reco_name << "\n";
                continue;
            }

            // Basic draw
            TCanvas c("c", "c", 800, 600);
            c.SetTicks(1,1);

            bool logy = LogAx(var, cfg);
            gPad->SetLogy(logy);

            // Style (keep minimal; you can customize)
            h_truth->SetLineWidth(2);
            h_truth->SetLineColor(kBlack);
            h_truth->SetStats(0);

            h_reco->SetLineWidth(2);
            h_reco->SetLineColor(kRed);
            h_reco->SetStats(0);

            h_truth->SetTitle((var + " truth vs reco").c_str());

            // Draw with sensible max
            const double maxy = std::max(h_truth->GetMaximum(), h_reco->GetMaximum());
            h_truth->SetMaximum(1.15 * maxy);

            h_truth->Draw("hist");
            h_reco->Draw("hist same");

            TLegend* leg = new TLegend(0.62, 0.7, 0.88, 0.85);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(h_truth, "Truth", "l");
            leg->AddEntry(h_reco,  "Reco",  "l");
            leg->Draw("same");

            const std::string out = outdir + var + "_truth_reco_compr" + wp_suffix + ".png";
            c.SaveAs(out.c_str());
        }
    }

    void PlotResponseMatrix() {
        if (!CheckFile()) return;

        const std::string outdir = GetDetRespOutDir();
        makeDirIfNeeded(outdir);

        for (const auto& var : pair_vars_for_det_resp) {
            const std::string h2_name =
                "h_" + var + "_vs_truth_" + var + "_single_b" + wp_filter;

            TH2* h2 = dynamic_cast<TH2*>(infile->Get(h2_name.c_str()));
            if (!h2) {
                std::cerr << "[PlotResponseMatrix] WARNING: missing " << h2_name << "\n";
                continue;
            }

            TCanvas c("c", "c", 850, 700);
            c.SetRightMargin(0.14);
            c.SetTicks(1,1);

            h2->SetTitle((var + " response matrix").c_str());
            h2->SetStats(0);
            h2->Draw("colz");

            const std::string out = outdir + var + "_response_matrix" + wp_suffix + ".png";
            c.SaveAs(out.c_str());
        }
    }

    void Plot1DRecoEffcy()
    {
        if (!CheckFile()) return;

        const std::string outdir = GetRecoEffcyOutDir();
        makeDirIfNeeded(outdir);

        for (const auto& var : pair_vars_for_1d_reco_effcy) {

            TCanvas c(Form("c_reco_eff1d_%s", var.c_str()), Form("c_reco_eff1d_%s", var.c_str()), 1400, 500);
            c.Divide(2,1);

            // ---- LEFT: same sign ----
            c.cd(1);
            gPad->SetTicks(1, 1);

            bool logx = LogAx(var, cfg);
            gPad->SetLogx(logx);

            const std::string hnum_name_ss   = "h_" + var + "_ss" + wp_filter;
            const std::string hdenom_name_ss = "h_" + var + "_ss" + truth_resn_filter;

            TH1D* h_num_ss   = dynamic_cast<TH1D*>(infile->Get(hnum_name_ss.c_str()));
            TH1D* h_den_ss   = dynamic_cast<TH1D*>(infile->Get(hdenom_name_ss.c_str()));

            if (!h_num_ss) {
                std::cerr << "[Plot1DRecoEffcy] WARNING: missing " << hnum_name_ss << "\n";
            }
            if (!h_den_ss) {
                std::cerr << "[Plot1DRecoEffcy] WARNING: missing " << hdenom_name_ss << "\n";
            }

            TGraphAsymmErrors* g_ss = nullptr;
            if (h_num_ss && h_den_ss) {
                g_ss = new TGraphAsymmErrors();
                g_ss->BayesDivide(h_num_ss, h_den_ss);

                g_ss->SetTitle("same sign");
                g_ss->GetXaxis()->SetTitle(h_den_ss->GetXaxis()->GetTitle());
                g_ss->GetYaxis()->SetTitle("#varepsilon");
                g_ss->GetYaxis()->SetRangeUser(0.0, 1.0);

                g_ss->Draw("AP");

                // // subplot title
                // TLatex lab;
                // lab.SetNDC();
                // lab.SetTextSize(0.045);
                // lab.DrawLatex(0.12, 0.92, "same sign");

                if (logx) adjustLogXRange(g_ss, h_den_ss);
            }

            // ---- RIGHT: opposite sign ----
            c.cd(2);
            gPad->SetTicks(1, 1);
            gPad->SetLogx(logx);

            const std::string hnum_name_op   = "h_" + var + "_op" + wp_filter;
            const std::string hdenom_name_op = "h_" + var + "_op" + truth_resn_filter;

            TH1D* h_num_op = dynamic_cast<TH1D*>(infile->Get(hnum_name_op.c_str()));
            TH1D* h_den_op = dynamic_cast<TH1D*>(infile->Get(hdenom_name_op.c_str()));

            if (!h_num_op) {
                std::cerr << "[Plot1DRecoEffcy] WARNING: missing " << hnum_name_op << "\n";
            }
            if (!h_den_op) {
                std::cerr << "[Plot1DRecoEffcy] WARNING: missing " << hdenom_name_op << "\n";
            }

            TGraphAsymmErrors* g_op = nullptr;
            if (h_num_op && h_den_op) {
                g_op = new TGraphAsymmErrors();
                g_op->BayesDivide(h_num_op, h_den_op);

                g_op->SetTitle("opposite sign");
                g_op->GetXaxis()->SetTitle(h_den_op->GetXaxis()->GetTitle());
                g_op->GetYaxis()->SetTitle("#varepsilon");
                g_op->GetYaxis()->SetRangeUser(0.0, 1.0);

                g_op->Draw("AP");

                TLatex lab;
                lab.SetNDC();
                lab.SetTextSize(0.045);
                lab.DrawLatex(0.12, 0.92, "opposite sign");

                if (logx) adjustLogXRange(g_op, h_den_op);
            }

            // save
            const std::string out = outdir + "reco_effcy_" + var + wp_suffix + ".png";
            c.SaveAs(out.c_str());
        }
    }

    void Plot2DRecoEffcySigned()
    {
        if (!CheckFile()) return;

        const std::string outdir = GetRecoEffcyOutDir();
        makeDirIfNeeded(outdir);

        for (const auto& var_pair : pair_vars_for_2d_reco_effcy) {
            const std::string& varx = var_pair[0];
            const std::string& vary = var_pair[1];

            TCanvas c(Form("c_reco_eff2d_%s_vs_%s", vary.c_str(), varx.c_str()), Form("c_reco_eff2d_%s_vs_%s", vary.c_str(), varx.c_str()), 1300, 550);
            c.Divide(2, 1);

            // "similar to Plot1DRecoEffcy": use logx decision via map; also allow logy for y-var if present
            const bool logx = LogAx(varx, cfg);
            const bool logy = LogAx(vary, cfg);

            // ---------- left pad: same sign ----------
            c.cd(1);
            gPad->SetTicks(1, 1);
            gPad->SetLogx(logx);
            gPad->SetLogy(logy);
            gPad->SetRightMargin(0.14);

            const std::string hnum_ss   = "h_" + vary + "_vs_" + varx + "_ss" + wp_filter;
            const std::string hden_ss   = "h_" + vary + "_vs_" + varx + "_ss" + truth_resn_filter;

            TH2D* h_num_ss = dynamic_cast<TH2D*>(infile->Get(hnum_ss.c_str()));
            TH2D* h_den_ss = dynamic_cast<TH2D*>(infile->Get(hden_ss.c_str()));

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

            const std::string hnum_op   = "h_" + vary + "_vs_" + varx + "_op" + wp_filter;
            const std::string hden_op   = "h_" + vary + "_vs_" + varx + "_op" + truth_resn_filter;

            TH2D* h_num_op = dynamic_cast<TH2D*>(infile->Get(hnum_op.c_str()));
            TH2D* h_den_op = dynamic_cast<TH2D*>(infile->Get(hden_op.c_str()));

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
        if (!CheckFile()) return;

        const std::string outdir = GetRecoEffcyOutDir();
        makeDirIfNeeded(outdir);

        for (const auto& var_pair : pair_vars_for_2d_reco_effcy) {
            const std::string& varx = var_pair[0];
            const std::string& vary = var_pair[1];

            TCanvas c(Form("c_reco_eff2d_%s_vs_%s_single_b", vary.c_str(), varx.c_str()),
                      Form("c_reco_eff2d_%s_vs_%s_single_b", vary.c_str(), varx.c_str()),
                      700, 600);
            c.cd();

            const bool logx = LogAx(varx, cfg);
            const bool logy = LogAx(vary, cfg);

            gPad->SetTicks(1, 1);
            gPad->SetLogx(logx);
            gPad->SetLogy(logy);
            gPad->SetRightMargin(0.14);

            const std::string hnum_single_b = "h_" + vary + "_vs_" + varx + "_single_b" + wp_filter;
            const std::string hden_single_b = "h_" + vary + "_vs_" + varx + "_single_b" + truth_resn_filter;

            TH2D* h_num_single_b = dynamic_cast<TH2D*>(infile->Get(hnum_single_b.c_str()));
            TH2D* h_den_single_b = dynamic_cast<TH2D*>(infile->Get(hden_single_b.c_str()));

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
    if (!CheckFile()) return;

    const std::string outdir = GetRecoEffcyOutDir();
    makeDirIfNeeded(outdir);

    for (const auto& var : pair_vars_for_1d_reco_effcy) {

        TCanvas c(Form("c_reco_eff1d_%s_single_b_op_compr", var.c_str()),
                  Form("c_reco_eff1d_%s_single_b_op_compr", var.c_str()),
                  800, 600);
        c.cd();
        gPad->SetTicks(1, 1);

        const bool logx = LogAx(var, cfg);
        gPad->SetLogx(logx);

        // ---- opposite sign ----
        const std::string hnum_name_op   = "h_" + var + "_op" + wp_filter;
        const std::string hdenom_name_op = "h_" + var + "_op" + truth_resn_filter;

        TH1D* h_num_op = dynamic_cast<TH1D*>(infile->Get(hnum_name_op.c_str()));
        TH1D* h_den_op = dynamic_cast<TH1D*>(infile->Get(hdenom_name_op.c_str()));

        if (!h_num_op)  { std::cerr << "[Plot1DRecoEffcySingleBOpCompr] WARNING: missing " << hnum_name_op   << "\n"; }
        if (!h_den_op)  { std::cerr << "[Plot1DRecoEffcySingleBOpCompr] WARNING: missing " << hdenom_name_op << "\n"; }

        // ---- single-b ----
        const std::string hnum_name_single_b   = "h_" + var + "_single_b" + wp_filter;
        const std::string hdenom_name_single_b = "h_" + var + "_single_b" + truth_resn_filter;

        TH1D* h_num_single_b = dynamic_cast<TH1D*>(infile->Get(hnum_name_single_b.c_str()));
        TH1D* h_den_single_b = dynamic_cast<TH1D*>(infile->Get(hdenom_name_single_b.c_str()));

        if (!h_num_single_b) { std::cerr << "[Plot1DRecoEffcySingleBOpCompr] WARNING: missing " << hnum_name_single_b   << "\n"; }
        if (!h_den_single_b) { std::cerr << "[Plot1DRecoEffcySingleBOpCompr] WARNING: missing " << hdenom_name_single_b << "\n"; }

        TGraphAsymmErrors* g_op = nullptr;
        TGraphAsymmErrors* g_single_b = nullptr;

        if (h_num_op && h_den_op) {
            g_op = new TGraphAsymmErrors();
            g_op->BayesDivide(h_num_op, h_den_op);

            g_op->GetXaxis()->SetTitle(h_den_op->GetXaxis()->GetTitle());
            g_op->SetLineColor(kBlack);
            g_op->SetMarkerColor(kBlack);
            g_op->SetLineWidth(2);
            g_op->SetMarkerStyle(20);
        }
        if (h_num_single_b && h_den_single_b) {
            g_single_b = new TGraphAsymmErrors();
            g_single_b->BayesDivide(h_num_single_b, h_den_single_b);

            g_single_b->GetXaxis()->SetTitle(h_den_single_b->GetXaxis()->GetTitle());
            g_single_b->SetLineColor(kRed);
            g_single_b->SetMarkerColor(kRed);
            g_single_b->SetLineWidth(2);
            g_single_b->SetMarkerStyle(20);
        }

        if (g_op || g_single_b) {
            // Prefer g_op as the frame if present; otherwise use g_single_b.
            TGraphAsymmErrors* g_frame = g_op ? g_op : g_single_b;
            g_frame->GetYaxis()->SetTitle("#varepsilon");
            g_frame->GetYaxis()->SetRangeUser(0.0, 1.0);
            g_frame->Draw("AP");

            if (g_op && g_op != g_frame)             g_op->Draw("P same");
            if (g_single_b && g_single_b != g_frame) g_single_b->Draw("P same");

            // Match log-x range behavior used in Plot1DRecoEffcy.
            // Use opposite-sign denom as the reference; if missing, fall back to single-b denom.
            if (logx) {
                if (h_den_op) adjustLogXRange(g_frame, h_den_op);
                else if (h_den_single_b) adjustLogXRange(g_frame, h_den_single_b);
            }

            TLegend* leg = new TLegend(0.62, 0.7, 0.85, 0.85);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            if (g_op)       leg->AddEntry(g_op,       "opposite sign", "lp");
            if (g_single_b) leg->AddEntry(g_single_b, "single-b",      "lp");
            leg->Draw("same");
        }

        const std::string out = outdir + "reco_effcy_" + var + "_single_b_op_compr" + wp_suffix + ".png";
        c.SaveAs(out.c_str());

        delete g_op;
        delete g_single_b;
    }
}

void Plot1DRecoEffcyRangedSingleBOpComprHelper(
    const std::string& varx,
    const std::string& vary,
    bool proj_axis,
    const std::vector<std::pair<float, float>>& proj_ranges)
{
    if (!CheckFile()) return;
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
    makeDirIfNeeded(outdir);
    const auto ranged_outdir = outdir + "ranged/";
    makeDirIfNeeded(ranged_outdir);

    const auto [nrows, ncols] = best_grid(static_cast<int>(proj_ranges.size()));
    TCanvas c(Form("c_reco_eff1d_ranged_%s_vs_%s%s", vary.c_str(), varx.c_str(), axisTag().c_str()),
              Form("c_reco_eff1d_ranged_%s_vs_%s%s", vary.c_str(), varx.c_str(), axisTag().c_str()),
              500 * ncols, 380 * nrows);
    c.Divide(ncols, nrows);

    const std::string xvar = proj_axis ? vary : varx;
    const bool logx = LogAx(xvar, cfg);

    for (std::size_t i = 0; i < proj_ranges.size(); ++i) {
        c.cd(static_cast<int>(i) + 1);
        gPad->SetTicks(1, 1);
        gPad->SetLogx(logx);

        const auto& range = proj_ranges[i];
        const std::string range_suffix = as_suffix(range);

        const std::string gname_op = "g_" + vary + "_vs_" + varx + "_op" + wp_filter + axisTag() + "_" + range_suffix + "_divided";
        const std::string gname_single_b = "g_" + vary + "_vs_" + varx + "_single_b" + wp_filter + axisTag() + "_" + range_suffix + "_divided";

        TGraphAsymmErrors* g_op = dynamic_cast<TGraphAsymmErrors*>(infile->Get(gname_op.c_str()));
        TGraphAsymmErrors* g_single_b = dynamic_cast<TGraphAsymmErrors*>(infile->Get(gname_single_b.c_str()));

        if (!g_op) {
            std::cerr << "[Plot1DRecoEffcyRangedSingleBOpComprHelper] WARNING: missing " << gname_op << "\n";
        }
        if (!g_single_b) {
            std::cerr << "[Plot1DRecoEffcyRangedSingleBOpComprHelper] WARNING: missing " << gname_single_b << "\n";
        }

        if (!g_op && !g_single_b) continue;

        if (g_op) {
            g_op->SetLineColor(kBlack);
            g_op->SetMarkerColor(kBlack);
            g_op->SetLineWidth(2);
            g_op->SetMarkerStyle(20);
        }
        if (g_single_b) {
            g_single_b->SetLineColor(kRed);
            g_single_b->SetMarkerColor(kRed);
            g_single_b->SetLineWidth(2);
            g_single_b->SetMarkerStyle(20);
        }

        TGraphAsymmErrors* g_frame = g_op ? g_op : g_single_b;
        g_frame->GetYaxis()->SetTitle("#varepsilon");
        g_frame->GetYaxis()->SetRangeUser(0.0, 1.0);
        g_frame->GetXaxis()->SetTitle(xvar.c_str());
        g_frame->Draw("AP");

        if (g_op && g_op != g_frame) g_op->Draw("P SAME");
        if (g_single_b && g_single_b != g_frame) g_single_b->Draw("P SAME");

        TLegend* leg = new TLegend(0.42, 0.62, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        if (g_op) leg->AddEntry(g_op, "opposite sign", "lp");
        if (g_single_b) leg->AddEntry(g_single_b, "single-b", "lp");

        const std::string ranged_var = proj_axis ? varx : vary;
        const std::string ranged_text = var_label(ranged_var) + ": [" + fmt_range_val(range.first) + ", " + fmt_range_val(range.second) + "]";
        leg->AddEntry((TObject*)nullptr, ranged_text.c_str(), "");
        leg->Draw("same");
    }

    const bool projx = !proj_axis;
    const std::string var_target = projx ? varx : vary;
    const std::string var_range = projx ? vary : varx;
    const std::string out = ranged_outdir + "reco_effcy_single_b_op_compr_" + var_target + "_in_" + var_range + "_bins" + wp_suffix + ".png";
    c.SaveAs(out.c_str());
}

void Plot1DRecoEffcyRangedSingleBOpCompr()
{
    if (!CheckFile()) return;

    for (const auto& cfg_tuple : reco_eff_proj_divide_cfgs) {
        const std::string& varx = std::get<0>(cfg_tuple);
        const std::string& vary = std::get<1>(cfg_tuple);
        const bool proj_axis = std::get<2>(cfg_tuple);
        const auto* proj_ranges = std::get<3>(cfg_tuple);

        if (proj_ranges == nullptr || proj_ranges->empty()) continue;
        Plot1DRecoEffcyRangedSingleBOpComprHelper(varx, vary, proj_axis, *proj_ranges);
    }
}

private:
    std::string GetDetRespOutDir() const {
        // Ensure trailing slash on data_dir (optional)
        std::string base = data_dir;
        if (!base.empty() && base.back() != '/') base += '/';
        return base + run_period + "_det_resp_plots" + no_reco_resn_cuts_file_dir_suffix + "/";
    }

    bool CheckFile() const {
        if (!infile || !infile->IsOpen()) {
            std::cerr << "[DetRespPlotter] ERROR: input file not open. Call Initialize().\n";
            return false;
        }
        return true;
    }

    std::string GetRecoEffcyOutDir() const
    {
        std::string base = data_dir;
        if (!base.empty() && base.back() != '/') base += '/';
        return base + run_period + "_reco_effcy_plots" + no_reco_resn_cuts_file_dir_suffix + "/";
    }

};
