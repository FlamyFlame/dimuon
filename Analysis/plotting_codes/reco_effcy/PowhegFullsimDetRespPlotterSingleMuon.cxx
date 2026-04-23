#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TError.h>
#include <TLatex.h>
#include <TPad.h>
#include <TGraphAsymmErrors.h>  // must precede TRatioPlot.h (which only forward-declares it)
#include <TRatioPlot.h>
#include <array>
#include <limits>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include "../../Utilities/PlotUtils.h"
#include "../../Utilities/PlotCommonConfig.h"

class PowhegFullsimDetRespPlotterSingleMuon {
public:

    PowhegFullsimDetRespPlotterSingleMuon(int run_year_input, bool tight_WP_input = false)
    :   run_year(run_year_input % 2000),
        tight_WP(tight_WP_input)
    {
        isRun3 = !(run_year <= 18);
    }

    PowhegFullsimDetRespPlotterSingleMuon()
    : PowhegFullsimDetRespPlotterSingleMuon(17, false) {}

    ~PowhegFullsimDetRespPlotterSingleMuon(){
        if (infile){
            infile->Close();
            delete infile;
            infile = nullptr;
        }
    }

    void Run(){
        if (!Initialize()) return;

        PlotTruthRecoCompr();
        PlotResponseMatrix();
        Plot1DRecoEffcy();
        Plot2DRecoEffcy();
        Plot1DRecoEffcySignCompr();
    }

protected:
    int  run_year = 17;
    bool isRun3;
    bool tight_WP;

    std::string run_period;
    std::string wp_suffix;
    std::string wp_filter;

    PlotCommonConfig cfg;

    std::string data_dir;
    std::string infile_name;
    TFile* infile{nullptr};

    // variables used for detector-response plots (reco vs truth comparisons)
    std::vector<std::string> muon_vars_for_det_resp = {"pt", "eta", "phi"};

    // variables used for 1-D reco-efficiency plots (truth-level x-axis)
    std::vector<std::string> muon_vars_for_1d_reco_effcy = {"truth_pt", "truth_eta", "truth_phi"};

    // 2-D reco-efficiency variable pairs (both truth-level)
    std::vector<std::array<std::string, 2>> muon_vars_for_2d_reco_effcy = {
        {"truth_pt",  "truth_eta"},
        {"truth_pt",  "truth_phi"},
        {"truth_eta", "truth_phi"}
    };

    // sign categories: _sign1 = mu+, _sign2 = mu-
    std::vector<std::string> signs      = {"_sign1", "_sign2"};
    std::vector<std::string> sign_labels = {"#mu^{+}", "#mu^{-}"};

    bool Initialize(){
        run_period = isRun3 ? "run3" : "run2";
        wp_suffix  = tight_WP ? "_tightWP" : "_mediumWP";
        wp_filter  = tight_WP ? "_pass_tight" : "_pass_medium";

        data_dir    = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample";
        infile_name = "histograms_powheg_fullsim_single_muon_pp" + std::to_string(run_year) + ".root";

        std::string infile_path = data_dir;
        if (!infile_path.empty() && infile_path.back() != '/') infile_path += '/';
        infile_path += infile_name;

        infile = TFile::Open(infile_path.c_str(), "READ");
        if (!infile || infile->IsZombie()){
            std::cerr << "[PowhegFullsimDetRespPlotterSingleMuon::Initialize] ERROR: failed to open file: "
                      << infile_path << "\n";
            if (infile){ delete infile; infile = nullptr; }
            return false;
        }
        return true;
    }

    // -----------------------------------------------------------------------
    // PlotTruthRecoCompr
    //   For each reco variable (pt / eta / phi) draw a 2-pad canvas with
    //   sign1 (left) and sign2 (right), each showing truth vs reco overlaid
    //   with a ratio panel (TRatioPlot).
    // -----------------------------------------------------------------------
    void PlotTruthRecoCompr(){
        if (!CheckFile()) return;

        const std::string outdir = GetDetRespOutDir();
        gSystem->mkdir(outdir.c_str(), kTRUE);

        for (const auto& var : muon_vars_for_det_resp){
            for (std::size_t is = 0; is < signs.size(); ++is){
                const std::string& sign = signs[is];

                const std::string h_truth_name = "h_truth_" + var + sign + wp_filter;
                const std::string h_reco_name  = "h_"       + var + sign + wp_filter;

                TH1* h_truth_in = GetHist<TH1>(h_truth_name);
                TH1* h_reco_in  = GetHist<TH1>(h_reco_name);

                if (!h_truth_in){
                    std::cerr << "[PlotTruthRecoCompr] WARNING: missing " << h_truth_name << "\n";
                    continue;
                }
                if (!h_reco_in){
                    std::cerr << "[PlotTruthRecoCompr] WARNING: missing " << h_reco_name << "\n";
                    continue;
                }

                TH1* h_truth = dynamic_cast<TH1*>(h_truth_in->Clone((h_truth_name + "_clone").c_str()));
                TH1* h_reco  = dynamic_cast<TH1*>(h_reco_in->Clone((h_reco_name  + "_clone").c_str()));
                if (!h_truth || !h_reco){
                    delete h_truth; delete h_reco; continue;
                }
                h_truth->SetDirectory(nullptr);
                h_reco->SetDirectory(nullptr);
                if (h_truth->GetSumw2N() == 0) h_truth->Sumw2();
                if (h_reco->GetSumw2N()  == 0) h_reco->Sumw2();

                TCanvas c("c", "c", 900, 760);
                c.cd();

                const bool logy = LogAx(var, cfg);

                h_truth->Scale(1, "width");
                h_reco ->Scale(1, "width");

                auto DiffY1D = [](const std::string& v) -> std::string {
                    if (v == "pt")  return "d#sigma/dp_{T} [pb/GeV]";
                    if (v == "eta") return "d#sigma/d#eta";
                    if (v == "phi") return "d#sigma/d#phi";
                    return "d#sigma/d(" + v + ")";
                };

                h_truth->SetLineWidth(2);
                h_truth->SetLineColor(kBlack);
                h_truth->SetStats(0);
                h_truth->SetTitle("");

                h_reco->SetLineWidth(2);
                h_reco->SetLineColor(kRed);
                h_reco->SetStats(0);
                h_reco->SetTitle("");

                TRatioPlot rp(h_reco, h_truth, "divsym");
                rp.SetLeftMargin(0.16);
                rp.SetRightMargin(0.05);
                rp.SetUpTopMargin(0.08);
                rp.SetUpBottomMargin(0.02);
                rp.SetLowTopMargin(0.02);
                rp.SetLowBottomMargin(0.40);
                rp.SetH1DrawOpt("hist");
                rp.SetH2DrawOpt("hist");
                rp.SetGraphDrawOpt("PE1");
                rp.Draw();
                rp.SetSplitFraction(0.30);
                rp.SetSeparationMargin(0.02);

                if (rp.GetUpperPad()){
                    rp.GetUpperPad()->SetTicks(1, 1);
                    rp.GetUpperPad()->SetLogy(logy);
                    rp.GetUpperPad()->cd();
                }

                if (rp.GetUpperRefYaxis()){
                    rp.GetUpperRefYaxis()->SetTitle(DiffY1D(var).c_str());
                    rp.GetUpperRefYaxis()->SetTitleSize(0.038);
                    rp.GetUpperRefYaxis()->SetLabelSize(0.030);
                    rp.GetUpperRefYaxis()->SetTitleOffset(1.30);
                }
                if (rp.GetUpperRefXaxis()){
                    rp.GetUpperRefXaxis()->SetTitleSize(0.0);
                    rp.GetUpperRefXaxis()->SetLabelSize(0.0);
                }
                if (rp.GetLowerRefYaxis()){
                    rp.GetLowerRefYaxis()->SetTitle("Reco/Truth");
                    rp.GetLowerRefYaxis()->SetNdivisions(4, kFALSE);
                    rp.GetLowerRefYaxis()->SetRangeUser(-1.2, 1.2);
                    rp.GetLowerRefYaxis()->SetTitleSize(0.030);
                    rp.GetLowerRefYaxis()->SetTitleOffset(1.20);
                    rp.GetLowerRefYaxis()->SetLabelSize(0.025);
                }
                if (rp.GetLowerRefXaxis()){
                    rp.GetLowerRefXaxis()->SetTitle(h_truth->GetXaxis()->GetTitle());
                    rp.GetLowerRefXaxis()->SetTitleSize(0.040);
                    rp.GetLowerRefXaxis()->SetTitleOffset(1.10);
                    rp.GetLowerRefXaxis()->SetLabelSize(0.032);
                }

                if (rp.GetUpperPad()) rp.GetUpperPad()->cd();

                TLegend* leg = new TLegend(0.62, 0.70, 0.88, 0.85);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->AddEntry(h_truth, "Truth", "l");
                leg->AddEntry(h_reco,  "Reco",  "l");
                leg->Draw("same");

                TLatex lab;
                lab.SetNDC();
                lab.SetTextSize(0.045);
                lab.DrawLatex(0.12, 0.92, (sign_labels[is] + " " + var + " truth vs reco").c_str());

                c.Modified(); c.Update();

                const std::string out = outdir + var + "_truth_reco_compr" + sign + wp_suffix + ".png";
                c.SaveAs(out.c_str());

                delete h_truth; delete h_reco;
            }
        }
    }

    // -----------------------------------------------------------------------
    // PlotResponseMatrix
    //   2-pad canvas (sign1 | sign2) for each {reco_var, truth_var} pair.
    // -----------------------------------------------------------------------
    void PlotResponseMatrix(){
        if (!CheckFile()) return;

        const std::string outdir = GetDetRespOutDir();
        gSystem->mkdir(outdir.c_str(), kTRUE);

        auto DiffZ2D = [](const std::string& v) -> std::string {
            if (v == "pt")  return "d^{2}#sigma/dp_{T}^{reco}dp_{T}^{truth} [pb/GeV^{2}]";
            if (v == "eta") return "d^{2}#sigma/d#eta^{reco}d#eta^{truth}";
            if (v == "phi") return "d^{2}#sigma/d#phi^{reco}d#phi^{truth}";
            return "d^{2}#sigma/d" + v + "^{reco}d" + v + "^{truth}";
        };

        for (const auto& var : muon_vars_for_det_resp){

            TCanvas c("c", "c", 1300, 550);
            c.Divide(2, 1);

            std::vector<TH2*> h2_clones;

            for (std::size_t is = 0; is < signs.size(); ++is){
                c.cd(static_cast<int>(is) + 1);
                gPad->SetRightMargin(0.14);
                gPad->SetTicks(1, 1);

                const std::string& sign = signs[is];
                const std::string h2_name = "h_" + var + "_vs_truth_" + var + sign + wp_filter;

                TH2* h2_in = GetHist<TH2>(h2_name);
                if (!h2_in){
                    std::cerr << "[PlotResponseMatrix] WARNING: missing " << h2_name << "\n";
                    continue;
                }
                TH2* h2 = dynamic_cast<TH2*>(h2_in->Clone(("h2_resp_" + var + sign).c_str()));
                h2->SetDirectory(nullptr);
                h2->Scale(1, "width");

                gPad->SetLogz(true);
                if (var == "pt"){
                    gPad->SetLogx(true);
                    gPad->SetLogy(true);
                }

                h2->SetTitle((sign_labels[is] + " " + var + " response matrix").c_str());
                h2->SetStats(0);
                h2->GetZaxis()->SetTitle(DiffZ2D(var).c_str());
                h2->Draw("colz");
                h2_clones.push_back(h2);
            }

            const std::string out = outdir + var + "_response_matrix" + wp_suffix + ".png";
            c.SaveAs(out.c_str());
            for (auto* h : h2_clones) delete h;
        }
    }

    // -----------------------------------------------------------------------
    // Plot1DRecoEffcy
    //   2-pad canvas (sign1 | sign2) for each truth variable.
    //   Efficiency = h_var_signX_wp_filter / h_var_signX
    // -----------------------------------------------------------------------
    void Plot1DRecoEffcy(){
        if (!CheckFile()) return;

        const std::string outdir = GetRecoEffcyOutDir();
        gSystem->mkdir(outdir.c_str(), kTRUE);

        for (const auto& var : muon_vars_for_1d_reco_effcy){

            TCanvas c(Form("c_reco_eff1d_%s", var.c_str()),
                      Form("c_reco_eff1d_%s", var.c_str()), 1400, 500);
            c.Divide(2, 1);

            const bool logx = LogAx(var, cfg);

            std::vector<TH1D*> h_effs; // collect clones; delete after SaveAs

            for (std::size_t is = 0; is < signs.size(); ++is){
                c.cd(static_cast<int>(is) + 1);
                gPad->SetTicks(1, 1);
                gPad->SetLogx(logx);

                const std::string& sign = signs[is];
                const std::string hnum_name   = "h_" + var + sign + wp_filter;
                const std::string hdenom_name = "h_" + var + sign;

                TH1D* h_num   = GetHist<TH1D>(hnum_name);
                TH1D* h_denom = GetHist<TH1D>(hdenom_name);

                if (!h_num)   std::cerr << "[Plot1DRecoEffcy] WARNING: missing " << hnum_name   << "\n";
                if (!h_denom) std::cerr << "[Plot1DRecoEffcy] WARNING: missing " << hdenom_name << "\n";

                if (h_num && h_denom){
                    TH1D* h_eff = dynamic_cast<TH1D*>(h_num->Clone(("h_eff_" + sign.substr(1) + "_" + var).c_str()));
                    h_eff->SetDirectory(nullptr);
                    h_eff->SetStats(0);
                    h_eff->Divide(h_denom);

                    h_eff->SetTitle(sign_labels[is].c_str());
                    h_eff->GetXaxis()->SetTitle(h_denom->GetXaxis()->GetTitle());
                    h_eff->GetYaxis()->SetTitle("#varepsilon");
                    h_eff->SetMinimum(0.0);
                    h_eff->SetMaximum(1.0);
                    h_eff->SetLineColor(kBlack);
                    h_eff->SetMarkerColor(kBlack);
                    h_eff->SetLineWidth(2);
                    h_eff->SetMarkerStyle(20);
                    h_eff->Draw("E1");

                    if (logx) AdjustLogXRangeForHist(h_eff, h_denom);
                    h_effs.push_back(h_eff);
                }
            }

            const std::string out = outdir + "reco_effcy_" + var + wp_suffix + ".png";
            c.SaveAs(out.c_str());
            for (auto* h : h_effs) delete h;
        }
    }

    // -----------------------------------------------------------------------
    // Plot2DRecoEffcy
    //   2-pad canvas (sign1 | sign2) for each 2-D truth-variable pair.
    // -----------------------------------------------------------------------
    void Plot2DRecoEffcy(){
        if (!CheckFile()) return;

        const std::string outdir = GetRecoEffcyOutDir();
        gSystem->mkdir(outdir.c_str(), kTRUE);

        for (const auto& var_pair : muon_vars_for_2d_reco_effcy){
            const std::string& varx = var_pair[0];
            const std::string& vary = var_pair[1];

            const bool logx = LogAx(varx, cfg);
            const bool logy = LogAx(vary, cfg);

            TCanvas c(Form("c_reco_eff2d_%s_vs_%s", vary.c_str(), varx.c_str()),
                      Form("c_reco_eff2d_%s_vs_%s", vary.c_str(), varx.c_str()), 1300, 550);
            c.Divide(2, 1);

            std::vector<TH2D*> h_effs; // collect clones; delete after SaveAs

            for (std::size_t is = 0; is < signs.size(); ++is){
                c.cd(static_cast<int>(is) + 1);
                gPad->SetTicks(1, 1);
                gPad->SetLogx(logx);
                gPad->SetLogy(logy);
                gPad->SetRightMargin(0.14);

                const std::string& sign  = signs[is];
                const std::string hnum_name   = "h_" + vary + "_vs_" + varx + sign + wp_filter;
                const std::string hdenom_name = "h_" + vary + "_vs_" + varx + sign;

                TH2D* h_num   = GetHist<TH2D>(hnum_name);
                TH2D* h_denom = GetHist<TH2D>(hdenom_name);

                if (!h_num)   std::cerr << "[Plot2DRecoEffcy] WARNING: missing " << hnum_name   << "\n";
                if (!h_denom) std::cerr << "[Plot2DRecoEffcy] WARNING: missing " << hdenom_name << "\n";

                if (h_num && h_denom){
                    TH2D* h_eff = dynamic_cast<TH2D*>(h_num->Clone(
                        ("h_eff_" + sign.substr(1) + "_" + vary + "_vs_" + varx).c_str()));
                    h_eff->SetDirectory(nullptr);
                    h_eff->SetStats(0);
                    h_eff->Divide(h_denom);
                    h_eff->SetMinimum(0.0);
                    h_eff->SetMaximum(1.0);
                    h_eff->SetTitle(sign_labels[is].c_str());
                    h_eff->Draw("colz");
                    h_effs.push_back(h_eff);
                }
            }

            const std::string out = outdir + "reco_effcy_" + vary + "_vs_" + varx + wp_suffix + ".png";
            c.SaveAs(out.c_str());
            for (auto* h : h_effs) delete h;
        }
    }

    // -----------------------------------------------------------------------
    // Plot1DRecoEffcySignCompr
    //   Single-pad canvas overlaying sign1 (mu+, black) and sign2 (mu-, red)
    //   efficiency curves for each truth variable.
    // -----------------------------------------------------------------------
    void Plot1DRecoEffcySignCompr(){
        if (!CheckFile()) return;

        const std::string outdir = GetRecoEffcyOutDir();
        gSystem->mkdir(outdir.c_str(), kTRUE);

        const std::vector<Color_t> colors = {kBlack, kRed};

        for (const auto& var : muon_vars_for_1d_reco_effcy){

            TCanvas c(Form("c_reco_eff1d_%s_sign_compr", var.c_str()),
                      Form("c_reco_eff1d_%s_sign_compr", var.c_str()), 800, 600);
            c.cd();
            gPad->SetTicks(1, 1);

            const bool logx = LogAx(var, cfg);
            gPad->SetLogx(logx);

            TLegend* leg = new TLegend(0.62, 0.70, 0.85, 0.85);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);

            TH1D* h_frame = nullptr;
            const TH1D* h_denom_ref = nullptr;
            std::vector<TH1D*> h_effs; // collect clones; delete after SaveAs

            for (std::size_t is = 0; is < signs.size(); ++is){
                const std::string& sign = signs[is];

                const std::string hnum_name   = "h_" + var + sign + wp_filter;
                const std::string hdenom_name = "h_" + var + sign;

                TH1D* h_num   = GetHist<TH1D>(hnum_name);
                TH1D* h_denom = GetHist<TH1D>(hdenom_name);

                if (!h_num)   std::cerr << "[Plot1DRecoEffcySignCompr] WARNING: missing " << hnum_name   << "\n";
                if (!h_denom) std::cerr << "[Plot1DRecoEffcySignCompr] WARNING: missing " << hdenom_name << "\n";

                if (!h_num || !h_denom) continue;

                TH1D* h_eff = dynamic_cast<TH1D*>(h_num->Clone(
                    ("h_eff_compr_" + sign.substr(1) + "_" + var).c_str()));
                h_eff->SetDirectory(nullptr);
                h_eff->SetStats(0);
                h_eff->Divide(h_denom);

                h_eff->GetXaxis()->SetTitle(h_denom->GetXaxis()->GetTitle());
                h_eff->SetLineColor(colors[is]);
                h_eff->SetMarkerColor(colors[is]);
                h_eff->SetLineWidth(2);
                h_eff->SetMarkerStyle(20);

                if (!h_frame){
                    h_eff->GetYaxis()->SetTitle("#varepsilon");
                    h_eff->SetMinimum(0.0);
                    h_eff->SetMaximum(1.0);
                    h_eff->SetTitle("");
                    h_eff->Draw("E1");
                    h_frame    = h_eff;
                    h_denom_ref = h_denom;
                } else {
                    h_eff->Draw("E1 same");
                }

                leg->AddEntry(h_eff, sign_labels[is].c_str(), "lp");
                h_effs.push_back(h_eff);
            }

            if (h_frame){
                if (logx && h_denom_ref) AdjustLogXRangeForHist(h_frame, h_denom_ref);
                leg->Draw("same");
            }

            const std::string out = outdir + "reco_effcy_" + var + "_sign_compr" + wp_suffix + ".png";
            c.SaveAs(out.c_str());
            for (auto* h : h_effs) delete h;
        }
    }

private:
    template <typename T>
    T* GetHist(const std::string& hname) const {
        if (!infile) return nullptr;
        return dynamic_cast<T*>(infile->Get(hname.c_str()));
    }

    bool CheckFile() const {
        if (!infile || !infile->IsOpen()){
            std::cerr << "[PowhegFullsimDetRespPlotterSingleMuon] ERROR: input file not open.\n";
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

    std::string GetDetRespOutDir() const {
        std::string base = data_dir;
        if (!base.empty() && base.back() != '/') base += '/';
        return base + "plots/" + run_period + "_det_resp_plots_single_muon/" + (tight_WP ? "tight/" : "medium/");
    }

    std::string GetRecoEffcyOutDir() const {
        std::string base = data_dir;
        if (!base.empty() && base.back() != '/') base += '/';
        return base + "plots/" + run_period + "_reco_effcy_plots_single_muon/" + (tight_WP ? "tight/" : "medium/");
    }
};
