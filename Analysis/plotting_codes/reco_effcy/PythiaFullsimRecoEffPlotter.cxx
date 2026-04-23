#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include <TRatioPlot.h>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <iostream>
#include "../../Utilities/PlotUtils.h"
#include "../../Utilities/PlotCommonConfig.h"
#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../MuonObjectsParamsAndHelpers/PbPbBaseClass.h"

class PythiaFullsimRecoEffPlotterBase {
public:

    PythiaFullsimRecoEffPlotterBase(bool tight_WP_input = false, bool require_signal_cuts_input = false)
    : tight_WP(tight_WP_input), require_signal_cuts(require_signal_cuts_input){}

    virtual ~PythiaFullsimRecoEffPlotterBase(){
        if (infile){ infile->Close(); delete infile; infile = nullptr; }
    }

    virtual void Run(){
        if (!Initialize()) return;
        RunForFilter("");
    }

    double min_denom_neff = 10.0;

protected:
    bool tight_WP;
    bool require_signal_cuts;

    std::string wp_suffix;
    std::string wp_filter;
    std::string require_signal_cuts_hist_suffix_num;
    std::string require_signal_cuts_hist_suffix_denom;
    std::string require_signal_cuts_file_dir_suffix;

    PlotCommonConfig cfg;

    TFile* infile{nullptr};

    virtual std::string GetDataDir() const = 0;
    virtual std::string GetPlotDirPrefix() const = 0;

    std::vector<std::string> pair_vars_for_det_resp = {
        "pair_pt", "minv_zoomin", "dr_zoomin"
    };
    std::vector<std::string> pair_vars_for_1d_reco_effcy = {
        "truth_pair_pt", "truth_pair_eta", "truth_dr_zoomin", "truth_dr_2_0", "truth_minv_zoomin"
    };
    std::vector<std::array<std::string,2>> pair_vars_for_2d_reco_effcy = {
        {"truth_pair_pt",   "truth_pair_eta"},
        {"truth_pair_pt",   "truth_dr_zoomin"},
        {"truth_pair_eta",  "truth_dr_zoomin"},
        {"truth_deta_zoomin", "truth_dphi_zoomin"},
        {"truth_pair_pt",   "truth_minv_zoomin"},
        {"truth_minv_zoomin", "truth_dr_zoomin"}
    };

    std::vector<std::pair<float,float>> dr_ranges_for_reco_effcy = {
        {0.0f, 0.2f}, {0.2f, 0.4f}, {0.4f, 0.6f}, {0.6f, 1.0f}
    };
    std::vector<std::pair<float,float>> pair_pT_ranges_for_reco_effcy_dR = {
        {8.0f, 12.0f}, {12.0f, 20.0f}, {20.0f, std::numeric_limits<float>::max()}
    };
    std::vector<std::tuple<std::string, std::string, bool, const std::vector<std::pair<float,float>>*>>
        reco_eff_proj_divide_cfgs = {
            {"truth_pair_pt",  "truth_dr_zoomin", false, &dr_ranges_for_reco_effcy},
            {"truth_pair_eta", "truth_dr_zoomin", false, &dr_ranges_for_reco_effcy},
            {"truth_pair_pt",  "truth_dr_zoomin", true,  &pair_pT_ranges_for_reco_effcy_dR}
        };

    void SuppressLowNeffBins(TH1D* h_eff, const TH1D* h_den) const {
        if (!h_eff || !h_den || min_denom_neff <= 0.0) return;
        for (int ib = 1; ib <= h_den->GetNbinsX(); ++ib){
            const double val = h_den->GetBinContent(ib);
            const double err = h_den->GetBinError(ib);
            if (err <= 0.0){
                if (val <= 0.0){ h_eff->SetBinContent(ib, 0.0); h_eff->SetBinError(ib, 0.0); }
                continue;
            }
            if ((val/err) * (val/err) < min_denom_neff){
                h_eff->SetBinContent(ib, 0.0);
                h_eff->SetBinError(ib, 0.0);
            }
        }
    }

    virtual std::string GetInputFilePath() const {
        return GetDataDir() + "histograms_pythia_fullsim_"
            + GetPlotDirPrefix() + "_no_data_resonance_cuts.root";
    }

    bool Initialize(){
        wp_suffix   = tight_WP ? "_tightWP"  : "_mediumWP";
        wp_filter   = tight_WP ? "_pass_tight" : "_pass_medium";
        require_signal_cuts_hist_suffix_num    = require_signal_cuts ? "_pass_signal_truth" : "";
        require_signal_cuts_hist_suffix_denom  = require_signal_cuts ? "_and_signal_truth_and_reco" : "";
        require_signal_cuts_file_dir_suffix    = require_signal_cuts ? "_require_signal_cuts" : "";

        const std::string infile_path = GetInputFilePath();

        infile = TFile::Open(infile_path.c_str(), "READ");
        if (!infile || infile->IsZombie()){
            std::cerr << "[PythiaFullsimRecoEffPlotterBase::Initialize] ERROR: failed to open "
                      << infile_path << "\n";
            if (infile){ delete infile; infile = nullptr; }
            return false;
        }
        return true;
    }

    // --- helpers ---

    std::string RecoEffNumFilter() const {
        if (!require_signal_cuts) return wp_filter;
        return wp_filter + require_signal_cuts_hist_suffix_denom;
    }
    std::string RecoEffDenomFilter() const {
        return require_signal_cuts ? require_signal_cuts_hist_suffix_num : "";
    }

    template <typename T>
    T* GetHist(const std::string& hname) const {
        if (!infile) return nullptr;
        return dynamic_cast<T*>(infile->Get(hname.c_str()));
    }

    std::string GetDetRespOutDir(const std::string& ctr_tag = "") const {
        return GetDataDir() + "plots/" + GetPlotDirPrefix() + "_det_resp_plots"
            + ctr_tag + "/" + (tight_WP ? "tight/" : "medium/");
    }
    std::string GetRecoEffcyOutDir(const std::string& ctr_tag = "") const {
        return GetDataDir() + "plots/" + GetPlotDirPrefix() + "_reco_effcy_plots"
            + require_signal_cuts_file_dir_suffix + ctr_tag + "/"
            + (tight_WP ? "tight/" : "medium/");
    }

    bool CheckFile() const {
        if (!infile || !infile->IsOpen()){
            std::cerr << "[PythiaFullsimRecoEffPlotterBase] ERROR: input file not open.\n";
            return false;
        }
        return true;
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

    void RunForFilter(const std::string& filter_prefix){
        const std::string ctr_tag = filter_prefix;
        PlotTruthRecoCompr(filter_prefix, ctr_tag);
        PlotResponseMatrix(filter_prefix, ctr_tag);
        Plot1DRecoEffcy(filter_prefix, ctr_tag);
        Plot2DRecoEffcySigned(filter_prefix, ctr_tag);
        Plot2DRecoEffcySingleB(filter_prefix, ctr_tag);
        Plot1DRecoEffcySingleBOpCompr(filter_prefix, ctr_tag);
        Plot1DRecoEffcyRangedSingleBOpCompr(filter_prefix, ctr_tag);
    }

    // --- plot functions ---
    // filter_prefix: inserted between pair_catgr and wp_filter in hist names
    //   PP: "" ; overlay: "_ctr0_5" etc.
    // ctr_tag: used for output directory differentiation

    void PlotTruthRecoCompr(const std::string& fp, const std::string& ctr_tag){
        if (!CheckFile()) return;
        const std::string outdir = GetDetRespOutDir(ctr_tag);
        gSystem->mkdir(outdir.c_str(), kTRUE);

        for (const auto& var : pair_vars_for_det_resp){
            const std::string h_truth_name = "h_truth_" + var + "_single_b" + fp + wp_filter;
            const std::string h_reco_name  = "h_"       + var + "_single_b" + fp + wp_filter;

            TH1* h_truth_in = GetHist<TH1>(h_truth_name);
            TH1* h_reco_in  = GetHist<TH1>(h_reco_name);

            if (!h_truth_in){ std::cerr << "[PlotTruthRecoCompr] WARNING: missing " << h_truth_name << "\n"; continue; }
            if (!h_reco_in) { std::cerr << "[PlotTruthRecoCompr] WARNING: missing " << h_reco_name  << "\n"; continue; }

            TH1* h_truth = dynamic_cast<TH1*>(h_truth_in->Clone((h_truth_name + "_clone").c_str()));
            TH1* h_reco  = dynamic_cast<TH1*>(h_reco_in ->Clone((h_reco_name  + "_clone").c_str()));
            if (!h_truth || !h_reco){ delete h_truth; delete h_reco; continue; }
            h_truth->SetDirectory(nullptr); h_reco->SetDirectory(nullptr);
            if (h_truth->GetSumw2N() == 0) h_truth->Sumw2();
            if (h_reco ->GetSumw2N() == 0) h_reco ->Sumw2();

            bool logy = LogAx(var, cfg);

            h_truth->Scale(1, "width");
            h_reco ->Scale(1, "width");

            auto DiffY1D = [](const std::string& v) -> std::string {
                if (v == "pair_pt")     return "d#sigma/dp_{T}^{pair} [pb/GeV]";
                if (v == "minv_zoomin") return "d#sigma/dm_{#mu#mu} [pb/GeV]";
                if (v == "dr_zoomin")   return "d#sigma/d#DeltaR";
                return "d#sigma/d(" + v + ")";
            };

            h_truth->SetLineWidth(2); h_truth->SetLineColor(kBlack); h_truth->SetStats(0); h_truth->SetTitle("");
            h_reco ->SetLineWidth(2); h_reco ->SetLineColor(kRed);   h_reco ->SetStats(0); h_reco ->SetTitle("");

            TCanvas c("c","c",900,760);
            TRatioPlot rp(h_reco, h_truth, "divsym");
            rp.SetLeftMargin(0.16); rp.SetRightMargin(0.05);
            rp.SetUpTopMargin(0.08); rp.SetUpBottomMargin(0.02);
            rp.SetLowTopMargin(0.02); rp.SetLowBottomMargin(0.40);
            rp.SetH1DrawOpt("hist"); rp.SetH2DrawOpt("hist");
            rp.SetGraphDrawOpt("PE1");
            rp.Draw();
            rp.SetSplitFraction(0.30); rp.SetSeparationMargin(0.02);

            if (rp.GetUpperPad()){ rp.GetUpperPad()->SetTicks(1,1); rp.GetUpperPad()->SetLogy(logy); rp.GetUpperPad()->cd(); }
            if (rp.GetUpperRefYaxis()){ rp.GetUpperRefYaxis()->SetTitle(DiffY1D(var).c_str()); rp.GetUpperRefYaxis()->SetTitleSize(0.038); rp.GetUpperRefYaxis()->SetLabelSize(0.030); rp.GetUpperRefYaxis()->SetTitleOffset(1.30); }
            if (rp.GetUpperRefXaxis()){ rp.GetUpperRefXaxis()->SetTitleSize(0.0); rp.GetUpperRefXaxis()->SetLabelSize(0.0); }
            if (rp.GetLowerRefYaxis()){ rp.GetLowerRefYaxis()->SetTitle("Reco/Truth"); rp.GetLowerRefYaxis()->SetNdivisions(4, kFALSE); rp.GetLowerRefYaxis()->SetRangeUser(-1.2, 1.2); rp.GetLowerRefYaxis()->SetTitleSize(0.030); rp.GetLowerRefYaxis()->SetTitleOffset(1.20); rp.GetLowerRefYaxis()->SetLabelSize(0.025); }
            if (rp.GetLowerRefXaxis()){ rp.GetLowerRefXaxis()->SetTitle(h_truth->GetXaxis()->GetTitle()); rp.GetLowerRefXaxis()->SetTitleSize(0.040); rp.GetLowerRefXaxis()->SetTitleOffset(1.10); rp.GetLowerRefXaxis()->SetLabelSize(0.032); }

            if (rp.GetUpperPad()) rp.GetUpperPad()->cd();
            TLegend* leg = new TLegend(0.62,0.7,0.88,0.85); leg->SetBorderSize(0); leg->SetFillStyle(0);
            leg->AddEntry(h_truth,"Truth","l"); leg->AddEntry(h_reco,"Reco","l"); leg->Draw("same");
            TLatex lab; lab.SetNDC(); lab.SetTextSize(0.045); lab.DrawLatex(0.12, 0.92, (var + " truth vs reco").c_str());

            c.Modified(); c.Update();
            c.SaveAs((outdir + var + "_truth_reco_compr" + wp_suffix + ".png").c_str());
            delete h_truth; delete h_reco;
        }
    }

    void PlotResponseMatrix(const std::string& fp, const std::string& ctr_tag){
        if (!CheckFile()) return;
        const std::string outdir = GetDetRespOutDir(ctr_tag);
        gSystem->mkdir(outdir.c_str(), kTRUE);

        for (const auto& var : pair_vars_for_det_resp){
            const std::string h2_name = "h_" + var + "_vs_truth_" + var + "_single_b" + fp + wp_filter;
            TH2* h2_in = GetHist<TH2>(h2_name);
            if (!h2_in){ std::cerr << "[PlotResponseMatrix] WARNING: missing " << h2_name << "\n"; continue; }
            TH2* h2 = dynamic_cast<TH2*>(h2_in->Clone(("h2_resp_" + var).c_str()));
            h2->SetDirectory(nullptr);
            h2->Scale(1, "width");

            auto DiffZ2D = [](const std::string& v) -> std::string {
                if (v == "pair_pt")     return "d^{2}#sigma/dp_{T}^{reco}dp_{T}^{truth} [pb/GeV^{2}]";
                if (v == "minv_zoomin") return "d^{2}#sigma/dm^{reco}dm^{truth} [pb/GeV^{2}]";
                if (v == "dr_zoomin")   return "d^{2}#sigma/d#DeltaR^{reco}d#DeltaR^{truth}";
                return "d^{2}#sigma/d" + v + "^{reco}d" + v + "^{truth}";
            };

            TCanvas c("c","c",850,700);
            c.SetRightMargin(0.14); c.SetTicks(1,1); c.cd();
            gPad->SetLogz(true);
            if (var == "pair_pt"){ gPad->SetLogx(true); gPad->SetLogy(true); }
            h2->SetStats(0); h2->SetTitle((var + " response matrix").c_str());
            h2->GetZaxis()->SetTitle(DiffZ2D(var).c_str());
            h2->Draw("colz");
            c.SaveAs((outdir + var + "_response_matrix" + wp_suffix + ".png").c_str());
            delete h2;
        }
    }

    void Plot1DRecoEffcy(const std::string& fp, const std::string& ctr_tag){
        if (!CheckFile()) return;
        const std::string outdir = GetRecoEffcyOutDir(ctr_tag);
        gSystem->mkdir(outdir.c_str(), kTRUE);

        const std::string num_filter   = RecoEffNumFilter();
        const std::string denom_filter = RecoEffDenomFilter();

        for (const auto& var : pair_vars_for_1d_reco_effcy){
            TCanvas c(Form("c_reco_eff1d_%s",var.c_str()), Form("c_reco_eff1d_%s",var.c_str()), 1400, 500);
            c.Divide(2,1);
            bool logx = LogAx(var, cfg);

            std::vector<TH1D*> effs_to_delete;
            for (int ipad = 0; ipad < 2; ++ipad){
                const std::string sign  = (ipad == 0) ? "_ss" : "_op";
                const std::string title = (ipad == 0) ? "same sign" : "opposite sign";

                c.cd(ipad+1);
                gPad->SetTicks(1,1); gPad->SetLogx(logx);

                const std::string hnum_name   = "h_" + var + sign + fp + num_filter;
                const std::string hdenom_name = "h_" + var + sign + fp + denom_filter;

                TH1D* h_num = GetHist<TH1D>(hnum_name);
                TH1D* h_den = GetHist<TH1D>(hdenom_name);

                if (!h_num){ std::cerr << "[Plot1DRecoEffcy] WARNING: missing " << hnum_name   << "\n"; }
                if (!h_den){ std::cerr << "[Plot1DRecoEffcy] WARNING: missing " << hdenom_name << "\n"; }

                if (h_num && h_den){
                    TH1D* h_eff = dynamic_cast<TH1D*>(h_num->Clone(("h_eff_" + sign.substr(1) + "_" + var).c_str()));
                    h_eff->SetDirectory(nullptr); h_eff->SetStats(0);
                    h_eff->Divide(h_den);
                    SuppressLowNeffBins(h_eff, h_den);
                    h_eff->SetTitle(title.c_str());
                    h_eff->GetXaxis()->SetTitle(h_den->GetXaxis()->GetTitle());
                    h_eff->GetYaxis()->SetTitle("#varepsilon");
                    h_eff->SetMinimum(0.0); h_eff->SetMaximum(1.0);
                    h_eff->SetLineColor(kBlack); h_eff->SetMarkerColor(kBlack);
                    h_eff->SetLineWidth(2); h_eff->SetMarkerStyle(20);
                    h_eff->Draw("E1");
                    if (logx) AdjustLogXRangeForHist(h_eff, h_den);
                    effs_to_delete.push_back(h_eff);
                }
            }

            c.SaveAs((outdir + "reco_effcy_" + var + "_ss_op" + wp_suffix + ".png").c_str());
            for (auto* h : effs_to_delete) delete h;
        }
    }

    void Plot2DRecoEffcySigned(const std::string& fp, const std::string& ctr_tag){
        if (!CheckFile()) return;
        const std::string outdir = GetRecoEffcyOutDir(ctr_tag);
        gSystem->mkdir(outdir.c_str(), kTRUE);

        const std::string num_filter   = RecoEffNumFilter();
        const std::string denom_filter = RecoEffDenomFilter();
        for (const auto& vars : pair_vars_for_2d_reco_effcy){
            const std::string& varx = vars[0];
            const std::string& vary = vars[1];

            const bool lx = LogAx(varx, cfg);
            const bool ly = LogAx(vary, cfg);

            TCanvas c(Form("c_reco_eff2d_%s_vs_%s", vary.c_str(), varx.c_str()),
                      Form("c_reco_eff2d_%s_vs_%s", vary.c_str(), varx.c_str()), 1300, 550);
            c.Divide(2, 1);

            TH2D* h_eff_ss = nullptr;
            TH2D* h_eff_op = nullptr;

            c.cd(1);
            gPad->SetTicks(1, 1); gPad->SetLogx(lx); gPad->SetLogy(ly); gPad->SetRightMargin(0.14);
            {
                const std::string hnum_name   = "h_" + vary + "_vs_" + varx + "_ss" + fp + num_filter;
                const std::string hdenom_name = "h_" + vary + "_vs_" + varx + "_ss" + fp + denom_filter;
                TH2D* h_num = GetHist<TH2D>(hnum_name);
                TH2D* h_den = GetHist<TH2D>(hdenom_name);
                if (!h_num) std::cerr << "[Plot2DRecoEffcySigned] WARNING: missing " << hnum_name   << "\n";
                if (!h_den) std::cerr << "[Plot2DRecoEffcySigned] WARNING: missing " << hdenom_name << "\n";
                if (h_num && h_den){
                    h_eff_ss = dynamic_cast<TH2D*>(h_num->Clone(("h_eff2d_ss_" + vary + "_vs_" + varx).c_str()));
                    h_eff_ss->SetDirectory(nullptr); h_eff_ss->SetStats(0);
                    h_eff_ss->Divide(h_den);
                    h_eff_ss->SetMinimum(0.0); h_eff_ss->SetMaximum(1.0);
                    h_eff_ss->SetTitle("same sign");
                    h_eff_ss->Draw("colz");
                }
            }

            c.cd(2);
            gPad->SetTicks(1, 1); gPad->SetLogx(lx); gPad->SetLogy(ly); gPad->SetRightMargin(0.14);
            {
                const std::string hnum_name   = "h_" + vary + "_vs_" + varx + "_op" + fp + num_filter;
                const std::string hdenom_name = "h_" + vary + "_vs_" + varx + "_op" + fp + denom_filter;
                TH2D* h_num = GetHist<TH2D>(hnum_name);
                TH2D* h_den = GetHist<TH2D>(hdenom_name);
                if (!h_num) std::cerr << "[Plot2DRecoEffcySigned] WARNING: missing " << hnum_name   << "\n";
                if (!h_den) std::cerr << "[Plot2DRecoEffcySigned] WARNING: missing " << hdenom_name << "\n";
                if (h_num && h_den){
                    h_eff_op = dynamic_cast<TH2D*>(h_num->Clone(("h_eff2d_op_" + vary + "_vs_" + varx).c_str()));
                    h_eff_op->SetDirectory(nullptr); h_eff_op->SetStats(0);
                    h_eff_op->Divide(h_den);
                    h_eff_op->SetMinimum(0.0); h_eff_op->SetMaximum(1.0);
                    h_eff_op->SetTitle("opposite sign");
                    h_eff_op->Draw("colz");
                }
            }

            c.SaveAs((outdir + "reco_effcy_2d_" + vary + "_vs_" + varx + wp_suffix + ".png").c_str());
            delete h_eff_ss; delete h_eff_op;
        }
    }

    void Plot2DRecoEffcySingleB(const std::string& fp, const std::string& ctr_tag){
        if (!CheckFile()) return;
        const std::string outdir = GetRecoEffcyOutDir(ctr_tag);
        gSystem->mkdir(outdir.c_str(), kTRUE);

        const std::string num_filter   = RecoEffNumFilter();
        const std::string denom_filter = RecoEffDenomFilter();

        for (const auto& vars : pair_vars_for_2d_reco_effcy){
            const std::string& varx = vars[0];
            const std::string& vary = vars[1];

            const std::string hnum_name   = "h_" + vary + "_vs_" + varx + "_single_b" + fp + num_filter;
            const std::string hdenom_name = "h_" + vary + "_vs_" + varx + "_single_b" + fp + denom_filter;

            TH2D* h_num = GetHist<TH2D>(hnum_name);
            TH2D* h_den = GetHist<TH2D>(hdenom_name);

            if (!h_num || !h_den){
                if (!h_num) std::cerr << "[Plot2DRecoEffcySingleB] WARNING: missing " << hnum_name   << "\n";
                if (!h_den) std::cerr << "[Plot2DRecoEffcySingleB] WARNING: missing " << hdenom_name << "\n";
                continue;
            }

            TH2D* h_eff = dynamic_cast<TH2D*>(h_num->Clone(("h_eff2d_" + vary + "_vs_" + varx + "_single_b").c_str()));
            h_eff->SetDirectory(nullptr); h_eff->SetStats(0);
            h_eff->Divide(h_den);
            h_eff->SetMinimum(0.0); h_eff->SetMaximum(1.0);

            TCanvas c("c","c",850,700);
            c.SetRightMargin(0.14); c.SetTicks(1,1); c.cd();
            if (LogAx(varx, cfg)) gPad->SetLogx(true);
            h_eff->SetTitle((vary + " vs " + varx + "_single_b reco effcy").c_str());
            h_eff->Draw("colz");
            c.SaveAs((outdir + "reco_effcy_2d_" + vary + "_vs_" + varx + "_single_b" + wp_suffix + ".png").c_str());
            delete h_eff;
        }
    }

    void Plot1DRecoEffcySingleBOpCompr(const std::string& fp, const std::string& ctr_tag){
        if (!CheckFile()) return;
        const std::string outdir = GetRecoEffcyOutDir(ctr_tag);
        gSystem->mkdir(outdir.c_str(), kTRUE);

        const std::string num_filter   = RecoEffNumFilter();
        const std::string denom_filter = RecoEffDenomFilter();

        for (const auto& var : pair_vars_for_1d_reco_effcy){
            TCanvas c(Form("c_reco_eff_compr_%s",var.c_str()), Form("c_reco_eff_compr_%s",var.c_str()), 700, 500);
            c.cd(); gPad->SetTicks(1,1);
            bool logx = LogAx(var, cfg); gPad->SetLogx(logx);

            TH1D* h_num_op       = GetHist<TH1D>("h_" + var + "_op"       + fp + num_filter);
            TH1D* h_den_op       = GetHist<TH1D>("h_" + var + "_op"       + fp + denom_filter);
            TH1D* h_num_single_b = GetHist<TH1D>("h_" + var + "_single_b" + fp + num_filter);
            TH1D* h_den_single_b = GetHist<TH1D>("h_" + var + "_single_b" + fp + denom_filter);

            TH1D* h_eff_op = nullptr, *h_eff_single_b = nullptr;
            if (h_num_op && h_den_op){
                h_eff_op = dynamic_cast<TH1D*>(h_num_op->Clone(("h_eff_op_" + var).c_str()));
                h_eff_op->SetDirectory(nullptr); h_eff_op->SetStats(0);
                h_eff_op->Divide(h_den_op); SuppressLowNeffBins(h_eff_op, h_den_op);
                h_eff_op->GetXaxis()->SetTitle(h_den_op->GetXaxis()->GetTitle());
                h_eff_op->SetLineColor(kBlack); h_eff_op->SetMarkerColor(kBlack);
                h_eff_op->SetLineWidth(2); h_eff_op->SetMarkerStyle(20);
            }
            if (h_num_single_b && h_den_single_b){
                h_eff_single_b = dynamic_cast<TH1D*>(h_num_single_b->Clone(("h_eff_single_b_" + var).c_str()));
                h_eff_single_b->SetDirectory(nullptr); h_eff_single_b->SetStats(0);
                h_eff_single_b->Divide(h_den_single_b); SuppressLowNeffBins(h_eff_single_b, h_den_single_b);
                h_eff_single_b->GetXaxis()->SetTitle(h_den_single_b->GetXaxis()->GetTitle());
                h_eff_single_b->SetLineColor(kRed); h_eff_single_b->SetMarkerColor(kRed);
                h_eff_single_b->SetLineWidth(2); h_eff_single_b->SetMarkerStyle(20);
            }

            if (h_eff_op || h_eff_single_b){
                TH1D* h_frame = h_eff_op ? h_eff_op : h_eff_single_b;
                h_frame->GetYaxis()->SetTitle("#varepsilon");
                h_frame->SetMinimum(0.0); h_frame->SetMaximum(1.0);
                h_frame->Draw("E1");
                if (h_eff_op       && h_eff_op       != h_frame) h_eff_op->Draw("E1 same");
                if (h_eff_single_b && h_eff_single_b != h_frame) h_eff_single_b->Draw("E1 same");
                if (logx){
                    if (h_den_op)       AdjustLogXRangeForHist(h_frame, h_den_op);
                    else if (h_den_single_b) AdjustLogXRangeForHist(h_frame, h_den_single_b);
                }
                TLegend* leg = new TLegend(0.62, 0.7, 0.85, 0.85);
                leg->SetBorderSize(0); leg->SetFillStyle(0);
                if (h_eff_op)       leg->AddEntry(h_eff_op,       "opposite sign", "lp");
                if (h_eff_single_b) leg->AddEntry(h_eff_single_b, "single-b",      "lp");
                leg->Draw("same");
            }

            c.SaveAs((outdir + "reco_effcy_" + var + "_single_b_op_compr" + wp_suffix + ".png").c_str());
            delete h_eff_op; delete h_eff_single_b;
        }
    }

    void Plot1DRecoEffcyRangedSingleBOpCompr(const std::string& fp, const std::string& ctr_tag){
        if (!CheckFile()) return;
        const std::string outdir = GetRecoEffcyOutDir(ctr_tag);
        const std::string ranged_outdir = outdir + "ranged/";
        gSystem->mkdir(ranged_outdir.c_str(), kTRUE);

        const std::string num_filter = RecoEffNumFilter();

        auto as_suffix = [](const std::pair<float,float>& range){
            auto fmt = [](float x){
                std::ostringstream oss; oss << std::fixed << std::setprecision(2) << x;
                std::string s = oss.str();
                for (auto& c : s) if (c == '.') c = '_';
                while (s.find('-') != std::string::npos) s.replace(s.find('-'),1,"minus");
                return s;
            };
            const bool upper_is_max = !std::isfinite(range.second)
                || range.second >= std::numeric_limits<float>::max() * 0.5f;
            if (upper_is_max) return fmt(range.first) + "_TO_MAX";
            return fmt(range.first) + "_TO_" + fmt(range.second);
        };

        auto fmt_range_val = [](float x){
            if (!std::isfinite(x) || x >= std::numeric_limits<float>::max() * 0.5f) return std::string("MAX");
            std::ostringstream oss; oss << std::fixed << std::setprecision(1) << x;
            std::string s = oss.str();
            if (s.find('.') != std::string::npos){
                while (!s.empty() && s.back() == '0') s.pop_back();
                if (!s.empty() && s.back() == '.') s.pop_back();
            }
            return s;
        };

        auto best_grid = [](int n) -> std::pair<int,int> {
            if (n <= 3) return {1, n};
            int br = n, bc = 1, bd = n-1, ba = n;
            for (int cols = 1; cols <= n; ++cols){
                int rows = (n + cols - 1) / cols;
                if (rows < cols) continue;
                int d = rows - cols, a = rows * cols;
                if (d < bd || (d == bd && a < ba)){ br = rows; bc = cols; bd = d; ba = a; }
            }
            return {br, bc};
        };

        for (const auto& proj_cfg : reco_eff_proj_divide_cfgs){
            const std::string& varx    = std::get<0>(proj_cfg);
            const std::string& vary    = std::get<1>(proj_cfg);
            const bool proj_axis       = std::get<2>(proj_cfg);
            const auto* proj_ranges    = std::get<3>(proj_cfg);
            if (!proj_ranges || proj_ranges->empty()) continue;

            const std::string axisTag = proj_axis ? "_py" : "_px";
            const std::string xvar    = proj_axis ? vary : varx;
            const bool logx           = LogAx(xvar, cfg);

            const auto [nrows, ncols] = best_grid(static_cast<int>(proj_ranges->size()));
            TCanvas c(Form("c_reco_eff1d_ranged_%s_vs_%s%s", vary.c_str(), varx.c_str(), axisTag.c_str()),
                      Form("c_reco_eff1d_ranged_%s_vs_%s%s", vary.c_str(), varx.c_str(), axisTag.c_str()),
                      500 * ncols, 380 * nrows);
            c.Divide(ncols, nrows);

            for (std::size_t i = 0; i < proj_ranges->size(); ++i){
                c.cd(static_cast<int>(i)+1);
                gPad->SetTicks(1,1); gPad->SetLogx(logx);

                const auto& range = (*proj_ranges)[i];
                const std::string range_suffix = as_suffix(range);
                const std::string label_range = fmt_range_val(range.first) + " - " + fmt_range_val(range.second);

                const std::string hname_op = "h_" + vary + "_vs_" + varx + "_op" + fp + num_filter + axisTag + "_" + range_suffix + "_divided";
                const std::string hname_sb = "h_" + vary + "_vs_" + varx + "_single_b" + fp + num_filter + axisTag + "_" + range_suffix + "_divided";

                TH1D* h_op = GetHist<TH1D>(hname_op);
                TH1D* h_sb = GetHist<TH1D>(hname_sb);
                if (!h_op) std::cerr << "[Plot1DRecoEffcyRangedSingleBOpCompr] WARNING: missing " << hname_op << "\n";
                if (!h_sb) std::cerr << "[Plot1DRecoEffcyRangedSingleBOpCompr] WARNING: missing " << hname_sb << "\n";
                if (!h_op && !h_sb) continue;

                if (h_op){ h_op->SetLineColor(kBlack); h_op->SetMarkerColor(kBlack); h_op->SetLineWidth(2); h_op->SetMarkerStyle(20); h_op->SetTitle(""); h_op->SetStats(0); }
                if (h_sb){ h_sb->SetLineColor(kRed);   h_sb->SetMarkerColor(kRed);   h_sb->SetLineWidth(2); h_sb->SetMarkerStyle(20); h_sb->SetTitle(""); h_sb->SetStats(0); }

                TH1D* h_frame = h_op ? h_op : h_sb;
                h_frame->GetYaxis()->SetTitle("#varepsilon");
                h_frame->SetMinimum(0.0); h_frame->SetMaximum(1.0);
                h_frame->Draw("E1");
                if (h_op && h_op != h_frame) h_op->Draw("E1 same");
                if (h_sb && h_sb != h_frame) h_sb->Draw("E1 same");
                if (logx) AdjustLogXRangeForHist(h_frame, h_frame);

                TLegend* leg = new TLegend(0.14, 0.75, 0.88, 0.88);
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
                if (h_op) leg->AddEntry(h_op, "opposite sign", "lp");
                if (h_sb) leg->AddEntry(h_sb, "single-b",      "lp");
                leg->Draw("same");

                TLatex lab; lab.SetNDC(); lab.SetTextSize(0.038);
                lab.DrawLatex(0.14, 0.70, label_range.c_str());
            }

            const std::string out = ranged_outdir + "reco_effcy_" + vary + "_vs_" + varx + "_single_b_op_compr"
                + axisTag + wp_suffix + ".png";
            c.SaveAs(out.c_str());
        }
    }

};


// ---------- PP child ----------
class PythiaFullsimRecoEffPlotter : public PythiaFullsimRecoEffPlotterBase {
public:
    PythiaFullsimRecoEffPlotter(bool tight_WP_input = false, bool require_signal_cuts_input = false)
    : PythiaFullsimRecoEffPlotterBase(tight_WP_input, require_signal_cuts_input){}
    ~PythiaFullsimRecoEffPlotter() override {}

protected:
    std::string GetDataDir() const override {
        return FullSimSampleInputDir(FullSimSampleType::pp);
    }
    std::string GetPlotDirPrefix() const override {
        return FullSimSamplePlotDir(FullSimSampleType::pp);
    }
};


// ---------- Overlay child ----------
class PythiaFullsimRecoEffPlotterOverlay
  : public PythiaFullsimRecoEffPlotterBase
  , public PbPbBaseClass<PythiaFullsimRecoEffPlotterOverlay>
{
    friend class PbPbBaseClass<PythiaFullsimRecoEffPlotterOverlay>;
    int RunYear() const { return 24; }

public:
    FullSimSampleType fullsim_sample_type;

    PythiaFullsimRecoEffPlotterOverlay(
        FullSimSampleType sample_type = FullSimSampleType::hijing,
        bool tight_WP_input = false,
        bool require_signal_cuts_input = false)
    : PythiaFullsimRecoEffPlotterBase(tight_WP_input, require_signal_cuts_input)
    , fullsim_sample_type(sample_type)
    {}
    ~PythiaFullsimRecoEffPlotterOverlay() override {}

    void Run() override {
        InitializePbPb();
        if (!Initialize()) return;

        RunForFilter("");

        for (int ictr = 0; ictr < nCtrBins; ++ictr){
            const std::string& ctr_tag = ctr_bins[ictr];
            RunForFilter(ctr_tag);
        }
    }

protected:
    std::string GetDataDir() const override {
        return FullSimSampleInputDir(fullsim_sample_type);
    }
    std::string GetPlotDirPrefix() const override {
        return FullSimSamplePlotDir(fullsim_sample_type);
    }
    std::string GetInputFilePath() const override {
        return GetDataDir() + "histograms_pythia_fullsim_"
            + FullSimSampleLabel(fullsim_sample_type) + "_no_data_resonance_cuts.root";
    }
};
