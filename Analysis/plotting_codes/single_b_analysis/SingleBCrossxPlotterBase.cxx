#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../helper_functions.c"
#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../../Utilities/proj_range_to_suffix.cxx"
#include "../../MuonObjectsParamsAndHelpers/DatasetTriggerMap.h"
#include "../../MuonObjectsParamsAndHelpers/PPBaseClass.h"

class SingleBCrossxPlotterBase {
protected:
    std::string input_file_path;
    std::string output_dir;
    TFile* fin{nullptr};

    CommonEffcyConfig cfg;
    std::vector<std::pair<float, float>> q_eta_bins;
    std::vector<std::pair<float, float>> dr_bins{{0.0f, 0.2f}, {0.2f, 0.4f}, {0.4f, 0.6f}, {0.6f, 1.0f}};

    std::array<int, 6> line_colors{{kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2}};
    std::array<int, 6> marker_styles{{20, 21, 22, 33, 34, 29}};

public:
    SingleBCrossxPlotterBase(const std::string& in_path, const std::string& out_dir)
        : input_file_path(in_path), output_dir(out_dir), q_eta_bins(cfg.q_eta_proj_ranges_coarse_incl_gap) {}

    virtual ~SingleBCrossxPlotterBase() {
        if (fin) {
            fin->Close();
            delete fin;
            fin = nullptr;
        }
    }

    bool Init() {
        // Check input file exists
        if (gSystem->AccessPathName(input_file_path.c_str())) {
            throw std::runtime_error("Input file does not exist: " + input_file_path);
        }

        gSystem->mkdir(output_dir.c_str(), true);
        fin = TFile::Open(input_file_path.c_str(), "READ");
        if (!fin || fin->IsZombie()) {
            throw std::runtime_error("Cannot open input file: " + input_file_path);
        }
        gStyle->SetOptStat(0);
        return true;
    }

    bool CheckHistogramExists(const std::string& hname, const std::string& htype = "TH1") {
        TObject* obj = GetHistObject(hname);
        if (!obj) {
            throw std::runtime_error("Missing histogram: " + hname);
        }
        if ((htype == "TH2D" && !dynamic_cast<TH2D*>(obj)) ||
            (htype == "TH3D" && !dynamic_cast<TH3D*>(obj))) {
            throw std::runtime_error("Histogram " + hname + " is not of type " + htype);
        }
        return true;
    }

    bool CheckHistogramNonEmpty(const std::string& hname) {
        TH1* h = dynamic_cast<TH1*>(GetHistObject(hname));
        if (!h) {
            throw std::runtime_error("Cannot retrieve histogram to check non-empty: " + hname);
        }
        if (h->GetEntries() == 0) {
            throw std::runtime_error("Histogram is empty (Entries=0): " + hname);
        }
        return true;
    }

    virtual void Run() = 0;

protected:
    virtual TObject* GetHistObject(const std::string& name) {
        if (!fin) throw std::runtime_error("Input file not opened, cannot retrieve: " + name);
        return fin->Get(name.c_str());
    }

    void Save2DColz(const std::string& hname, const std::string& png_name,
                    const std::string& z_title = "") {
        CheckHistogramExists(hname, "TH2D");
        CheckHistogramNonEmpty(hname);

        TH2D* h = dynamic_cast<TH2D*>(GetHistObject(hname));
        if (!h) {
            throw std::runtime_error("Failed to retrieve TH2D: " + hname);
        }
        // Clone so we don't modify the in-file object
        TH2D* hplot = dynamic_cast<TH2D*>(h->Clone((hname + "_plot").c_str()));
        hplot->SetDirectory(nullptr);
        // Patch axis titles to physics notation
        {
            auto fix = [](TAxis* ax) {
                const std::string t = ax->GetTitle();
                if (t == "p_{T} [GeV]")  ax->SetTitle("p_{T}^{pair} [GeV]");
                else if (t == "#eta")     ax->SetTitle("#eta^{pair}");
            };
            fix(hplot->GetXaxis());
            fix(hplot->GetYaxis());
        }
        // Make differential: divide by bin widths of both axes
        hplot->Scale(1., "width");

        TCanvas c("c2d", "c2d", 900, 750);
        c.cd();
        c.SetRightMargin(0.17);
        c.SetLeftMargin(0.12);
        c.SetBottomMargin(0.12);

        const std::string xt = hplot->GetXaxis()->GetTitle();
        const std::string yt = hplot->GetYaxis()->GetTitle();
        const bool x_is_pair_pt = xt.find("p_{T}") != std::string::npos;
        const bool y_is_pair_pt = yt.find("p_{T}") != std::string::npos;
        if (x_is_pair_pt) c.SetLogx();
        if (y_is_pair_pt) c.SetLogy();
        if (x_is_pair_pt || y_is_pair_pt) c.SetLogz();

        hplot->GetXaxis()->SetTitleSize(0.05);
        hplot->GetYaxis()->SetTitleSize(0.05);
        hplot->GetZaxis()->SetTitleSize(0.045);
        hplot->GetXaxis()->SetLabelSize(0.04);
        hplot->GetYaxis()->SetLabelSize(0.04);
        hplot->GetZaxis()->SetLabelSize(0.035);
        hplot->GetYaxis()->SetTitleOffset(1.35);

        hplot->GetZaxis()->SetTitle(z_title.empty() ? "d^{2}N/d[X]d[Y]" : z_title.c_str());
        hplot->GetZaxis()->SetTitleOffset(1.35);
        hplot->Draw("colz");
        std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;
        delete hplot;
    }

    void DrawPairPtByEtaWithDrLines(
        const std::string& h3_name,
        const std::string& data_info_line1,
        const std::string& data_info_line2,
        const std::string& png_name,
        const std::string& y_title = "d#sigma/dp_{T} [pb GeV^{-1}]")
    {
        CheckHistogramExists(h3_name, "TH3D");
        CheckHistogramNonEmpty(h3_name);

        TH3D* h3 = dynamic_cast<TH3D*>(GetHistObject(h3_name));
        if (!h3) {
            throw std::runtime_error("Failed to retrieve TH3D: " + h3_name);
        }

        int nrow = 1;
        int ncol = 1;
        DetermineSubplotGrid(static_cast<int>(q_eta_bins.size()), nrow, ncol);

        TCanvas c("cpt", "pair_pt by eta with dr lines", 450 * ncol, 350 * nrow);
        c.Divide(ncol, nrow);

        std::vector<std::vector<TH1D*>> all_lines(q_eta_bins.size());
        std::vector<TLegend*> all_legends;
        all_legends.reserve(q_eta_bins.size());

        for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
            c.cd(static_cast<int>(ieta) + 1);
            gPad->SetLogx();
            gPad->SetLogy();
            gPad->SetLeftMargin(0.16);
            gPad->SetBottomMargin(0.13);

            const auto& eta_bin = q_eta_bins.at(ieta);
            const int y1 = h3->GetYaxis()->FindBin(eta_bin.first + 1e-6);
            const int y2 = h3->GetYaxis()->FindBin(eta_bin.second - 1e-6);

            double max_y = 0.0;
            std::vector<TH1D*> lines;
            lines.reserve(dr_bins.size());

            for (size_t idr = 0; idr < dr_bins.size(); ++idr) {
                const auto& dr_bin = dr_bins.at(idr);
                const int z1 = h3->GetZaxis()->FindBin(dr_bin.first + 1e-6);
                const int z2 = h3->GetZaxis()->FindBin(dr_bin.second - 1e-6);

                const std::string proj_name = "hpt_eta" + std::to_string(ieta) + "_dr" + std::to_string(idr) + "_" + std::to_string(std::rand());
                TH1D* hp = h3->ProjectionX(proj_name.c_str(), y1, y2, z1, z2, "e");
                hp->SetDirectory(nullptr);
                hp->Scale(1.0, "width");
                hp->SetLineWidth(2);
                hp->SetLineColor(line_colors.at(idr % line_colors.size()));
                hp->SetMarkerColor(line_colors.at(idr % line_colors.size()));
                hp->SetMarkerStyle(marker_styles.at(idr % marker_styles.size()));
                hp->SetMarkerSize(0.9);
                hp->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
                hp->GetYaxis()->SetTitle(y_title.c_str());
                hp->GetXaxis()->SetTitleSize(0.06);
                hp->GetYaxis()->SetTitleSize(0.06);
                hp->GetXaxis()->SetLabelSize(0.05);
                hp->GetYaxis()->SetLabelSize(0.05);
                hp->GetYaxis()->SetTitleOffset(1.45);
                hp->SetTitle("");

                max_y = std::max(max_y, hp->GetMaximum());
                lines.push_back(hp);
            }

            if (lines.empty()) continue;

            lines.at(0)->SetMaximum(max_y * 1.45);
            lines.at(0)->Draw("E1");
            for (size_t il = 1; il < lines.size(); ++il) {
                lines.at(il)->Draw("E1 SAME");
            }

            // Info text (no symbols): right-aligned, minimal margin
            TLegend* leg_info = new TLegend(0.38, 0.72, 0.93, 0.90);
            leg_info->SetBorderSize(0);
            leg_info->SetFillStyle(0);
            leg_info->SetTextSize(0.045);
            leg_info->SetTextAlign(32);
            leg_info->SetMargin(0.01);
            leg_info->AddEntry((TObject*)0, data_info_line1.c_str(), "");
            leg_info->AddEntry((TObject*)0, data_info_line2.c_str(), "");
            leg_info->AddEntry((TObject*)0, Form("#eta^{pair} #in [%.1f, %.1f]", eta_bin.first, eta_bin.second), "");
            leg_info->Draw();
            all_legends.push_back(leg_info);

            // dR lines: narrow legend, left-aligned text, symbol just left of text
            TLegend* leg_dr = new TLegend(0.70, 0.48, 0.93, 0.72);
            leg_dr->SetBorderSize(0);
            leg_dr->SetFillStyle(0);
            leg_dr->SetTextSize(0.045);
            leg_dr->SetTextAlign(12);
            leg_dr->SetMargin(0.22);
            for (size_t idr = 0; idr < dr_bins.size(); ++idr) {
                const auto& dr_bin = dr_bins.at(idr);
                leg_dr->AddEntry(lines.at(idr), Form("#DeltaR #in [%.1f, %.1f]", dr_bin.first, dr_bin.second), "lep");
            }
            leg_dr->Draw();
            all_legends.push_back(leg_dr);

            all_lines.at(ieta) = std::move(lines);
        }

        std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;

        for (auto& lines : all_lines) {
            for (TH1D* h : lines) delete h;
        }
        for (TLegend* leg : all_legends) delete leg;
    }

    // dR-integrated pair-pT vs pair-eta: one subplot per eta bin, single line per subplot.
    // h2_name: TH2D with pair_pt on X, pair_eta on Y.
    // differential=true: divide by bin width (d/dp_T); false: raw counts per eta bin.
    void DrawPairPtByEta(
        const std::string& h2_name,
        const std::string& data_info_line1,
        const std::string& data_info_line2,
        const std::string& png_name,
        const std::string& y_title = "d#sigma/dp_{T} [pb GeV^{-1}]",
        bool differential = true)
    {
        CheckHistogramExists(h2_name, "TH2D");
        CheckHistogramNonEmpty(h2_name);

        TH2D* h2 = dynamic_cast<TH2D*>(GetHistObject(h2_name));
        if (!h2) {
            throw std::runtime_error("Failed to retrieve TH2D: " + h2_name);
        }

        int nrow = 1, ncol = 1;
        DetermineSubplotGrid(static_cast<int>(q_eta_bins.size()), nrow, ncol);

        TCanvas c("cpt_eta", "pair_pt by eta dR-integrated", 450 * ncol, 350 * nrow);
        c.Divide(ncol, nrow);

        std::vector<TH1D*> all_hists;
        all_hists.reserve(q_eta_bins.size());
        std::vector<TLegend*> all_legends;
        all_legends.reserve(q_eta_bins.size());

        for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
            c.cd(static_cast<int>(ieta) + 1);
            gPad->SetLogx();
            gPad->SetLogy();
            gPad->SetLeftMargin(0.16);
            gPad->SetBottomMargin(0.13);

            const auto& eta_bin = q_eta_bins.at(ieta);
            const int y1 = h2->GetYaxis()->FindBin(eta_bin.first  + 1e-6);
            const int y2 = h2->GetYaxis()->FindBin(eta_bin.second - 1e-6);

            const std::string proj_name = "hpt_eta" + std::to_string(ieta) + "_"
                                          + std::to_string(std::rand());
            TH1D* hp = h2->ProjectionX(proj_name.c_str(), y1, y2, "e");
            hp->SetDirectory(nullptr);
            if (differential) hp->Scale(1.0, "width");
            hp->SetLineWidth(2);
            hp->SetLineColor(kBlack);
            hp->SetMarkerColor(kBlack);
            hp->SetMarkerStyle(20);
            hp->SetMarkerSize(0.9);
            hp->GetXaxis()->SetTitle("p_{T}^{pair} [GeV]");
            hp->GetYaxis()->SetTitle(y_title.c_str());
            hp->GetXaxis()->SetTitleSize(0.06);
            hp->GetYaxis()->SetTitleSize(0.06);
            hp->GetXaxis()->SetLabelSize(0.05);
            hp->GetYaxis()->SetLabelSize(0.05);
            hp->GetYaxis()->SetTitleOffset(1.45);
            hp->SetTitle("");

            if (hp->GetMaximum() > 0) hp->SetMaximum(hp->GetMaximum() * 1.45);
            hp->Draw("E1");

            TLegend* leg = new TLegend(0.60, 0.70, 0.93, 0.90);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.045);
            leg->SetTextAlign(32);
            leg->SetMargin(0.06);
            leg->AddEntry((TObject*)0, data_info_line1.c_str(), "");
            leg->AddEntry((TObject*)0, data_info_line2.c_str(), "");
            leg->AddEntry((TObject*)0, Form("#eta^{pair} #in [%.1f, %.1f]", eta_bin.first, eta_bin.second), "");
            leg->Draw();
            all_legends.push_back(leg);
            all_hists.push_back(hp);
        }

        std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;

        for (TH1D* h : all_hists) delete h;
        for (TLegend* leg : all_legends) delete leg;
    }
};
