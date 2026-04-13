#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../helper_functions.c"
#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"

class SignalAcceptancePlotter {
protected:
    std::string input_file_path;
    std::string output_dir;
    TFile* fin{nullptr};

    CommonEffcyConfig cfg;
    std::vector<std::pair<float, float>> q_eta_bins;

public:
    bool use_pt_bins_150 = false;

    SignalAcceptancePlotter(const std::string& in_path, const std::string& out_dir)
        : input_file_path(in_path), output_dir(out_dir),
          q_eta_bins(cfg.q_eta_proj_ranges_coarse_incl_gap) {}

    virtual ~SignalAcceptancePlotter() {
        if (fin) { fin->Close(); delete fin; fin = nullptr; }
    }

    bool Init() {
        if (gSystem->AccessPathName(input_file_path.c_str())) {
            throw std::runtime_error("Input file not found: " + input_file_path);
        }
        gSystem->mkdir(output_dir.c_str(), true);
        fin = TFile::Open(input_file_path.c_str(), "READ");
        if (!fin || fin->IsZombie()) {
            throw std::runtime_error("Cannot open input file: " + input_file_path);
        }
        gStyle->SetOptStat(0);
        return true;
    }

    // 2D colz plot of the acceptance ratio histogram.
    void Draw2DColz(const std::string& png_name) {
        const std::string hname = use_pt_bins_150 ? "h2d_sig_accept_pt_150_eta" : "h2d_sig_accept_pt_eta";
        TH2D* h = dynamic_cast<TH2D*>(fin->Get(hname.c_str()));
        if (!h) throw std::runtime_error("Missing histogram: " + hname);
        if (h->GetEntries() == 0) throw std::runtime_error(hname + " is empty");

        TH2D* hplot = dynamic_cast<TH2D*>(h->Clone((hname + "_plot").c_str()));
        hplot->SetDirectory(nullptr);

        TCanvas c("c2d", "", 900, 750);
        c.SetRightMargin(0.17);
        c.SetLeftMargin(0.12);
        c.SetBottomMargin(0.12);
        c.SetLogx();

        hplot->GetXaxis()->SetTitleSize(0.05);
        hplot->GetYaxis()->SetTitleSize(0.05);
        hplot->GetZaxis()->SetTitleSize(0.045);
        hplot->GetXaxis()->SetLabelSize(0.04);
        hplot->GetYaxis()->SetLabelSize(0.04);
        hplot->GetZaxis()->SetLabelSize(0.035);
        hplot->GetYaxis()->SetTitleOffset(1.35);
        hplot->GetZaxis()->SetTitleOffset(1.35);
        hplot->GetZaxis()->SetTitle("#alpha_{signal}");
        hplot->SetMinimum(0.);
        hplot->SetMaximum(1.);
        hplot->Draw("colz");

        const std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;
        delete hplot;
    }

    // 1D acceptance vs truth pair pT, one subplot per coarse eta bin.
    // Projects num and denom 2D histograms per eta bin, divides bin-by-bin.
    void DrawAcceptancePtByEta(const std::string& png_name,
                                const std::string& mc_label_line1,
                                const std::string& mc_label_line2) {
        const std::string num_name = use_pt_bins_150 ? "h2d_sig_accept_num_pt_150_eta" : "h2d_sig_accept_num_pt_eta";
        const std::string den_name = use_pt_bins_150 ? "h2d_sig_accept_denom_pt_150_eta" : "h2d_sig_accept_denom_pt_eta";
        TH2D* hnum = dynamic_cast<TH2D*>(fin->Get(num_name.c_str()));
        TH2D* hden = dynamic_cast<TH2D*>(fin->Get(den_name.c_str()));
        if (!hnum) throw std::runtime_error("Missing histogram: " + num_name);
        if (!hden) throw std::runtime_error("Missing histogram: " + den_name);
        if (hden->GetEntries() == 0) throw std::runtime_error(den_name + " is empty");

        int nrow = 1, ncol = 1;
        DetermineSubplotGrid(static_cast<int>(q_eta_bins.size()), nrow, ncol);

        TCanvas c("cacc", "signal acceptance by eta", 450 * ncol, 350 * nrow);
        c.Divide(ncol, nrow);

        std::vector<TH1D*> all_ratios;
        std::vector<TLegend*> all_legs;
        all_legs.reserve(q_eta_bins.size());

        for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
            c.cd(static_cast<int>(ieta) + 1);
            gPad->SetLogx();
            gPad->SetLeftMargin(0.16);
            gPad->SetBottomMargin(0.13);

            const auto& eta_bin = q_eta_bins.at(ieta);
            const int y1 = hnum->GetYaxis()->FindBin(eta_bin.first  + 1e-6f);
            const int y2 = hnum->GetYaxis()->FindBin(eta_bin.second - 1e-6f);

            const std::string tag = std::to_string(ieta);
            TH1D* hn = hnum->ProjectionX(("hn_acc_eta" + tag).c_str(), y1, y2, "e");
            TH1D* hd = hden->ProjectionX(("hd_acc_eta" + tag).c_str(), y1, y2, "e");
            hn->SetDirectory(nullptr);
            hd->SetDirectory(nullptr);

            TH1D* hratio = static_cast<TH1D*>(hn->Clone(("hratio_acc_eta" + tag).c_str()));
            hratio->SetDirectory(nullptr);
            hratio->Divide(hd);

            hratio->SetLineWidth(2);
            hratio->SetLineColor(kBlack);
            hratio->SetMarkerColor(kBlack);
            hratio->SetMarkerStyle(20);
            hratio->SetMarkerSize(0.9);
            hratio->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
            hratio->GetYaxis()->SetTitle("#alpha_{signal}");
            hratio->GetXaxis()->SetTitleSize(0.06);
            hratio->GetYaxis()->SetTitleSize(0.06);
            hratio->GetXaxis()->SetLabelSize(0.05);
            hratio->GetYaxis()->SetLabelSize(0.05);
            hratio->GetYaxis()->SetTitleOffset(1.45);
            hratio->SetTitle("");
            hratio->SetMinimum(0.);
            hratio->SetMaximum(1.3);
            hratio->Draw("E1");

            TLegend* leg = new TLegend(0.18, 0.65, 0.62, 0.92);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.045);
            leg->SetTextAlign(12);
            leg->SetMargin(0.06);
            leg->AddEntry((TObject*)0, mc_label_line1.c_str(), "");
            leg->AddEntry((TObject*)0, mc_label_line2.c_str(), "");
            leg->AddEntry((TObject*)0, Form("#eta^{pair} #in [%.1f, %.1f]", eta_bin.first, eta_bin.second), "");
            leg->Draw();
            all_legs.push_back(leg);

            all_ratios.push_back(hratio);
            delete hn;
            delete hd;
        }

        const std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;

        for (TH1D* h : all_ratios) delete h;
        for (TLegend* l : all_legs) delete l;
    }

    virtual void Run() = 0;
};
