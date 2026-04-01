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
        if (!fin) {
            throw std::runtime_error("Input file not opened");
        }
        TObject* obj = fin->Get(hname.c_str());
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
        TH1* h = dynamic_cast<TH1*>(fin->Get(hname.c_str()));
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
    void Save2DColz(const std::string& hname, const std::string& png_name) {
        CheckHistogramExists(hname, "TH2D");
        CheckHistogramNonEmpty(hname);

        TH2D* h = dynamic_cast<TH2D*>(fin->Get(hname.c_str()));
        if (!h) {
            throw std::runtime_error("Failed to retrieve TH2D: " + hname);
        }

        TCanvas c("c2d", "c2d", 800, 700);
        c.cd();
        h->Draw("colz");
        std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;
    }

    void DrawPairPtByEtaWithDrLines(
        const std::string& h3_name,
        const std::string& data_spec,
        const std::string& png_name)
    {
        CheckHistogramExists(h3_name, "TH3D");
        CheckHistogramNonEmpty(h3_name);

        TH3D* h3 = dynamic_cast<TH3D*>(fin->Get(h3_name.c_str()));
        if (!h3) {
            throw std::runtime_error("Failed to retrieve TH3D: " + h3_name);
        }

        int nrow = 1;
        int ncol = 1;
        DetermineSubplotGrid(static_cast<int>(q_eta_bins.size()), nrow, ncol);

        TCanvas c("cpt", "pair_pt by eta with dr lines", 450 * ncol, 350 * nrow);
        c.Divide(ncol, nrow);

        for (size_t ieta = 0; ieta < q_eta_bins.size(); ++ieta) {
            c.cd(static_cast<int>(ieta) + 1);

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
                hp->SetLineWidth(2);
                hp->SetLineColor(line_colors.at(idr % line_colors.size()));
                hp->SetMarkerColor(line_colors.at(idr % line_colors.size()));
                hp->SetMarkerStyle(marker_styles.at(idr % marker_styles.size()));
                hp->SetMarkerSize(0.9);
                hp->GetXaxis()->SetTitle("pair p_{T} [GeV]");
                hp->GetYaxis()->SetTitle("d#sigma/dp_{T}");
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

            TLegend leg(0.43, 0.52, 0.90, 0.90);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.AddEntry((TObject*)0, data_spec.c_str(), "");
            leg.AddEntry((TObject*)0, Form("q#eta #in [%.1f, %.1f]", eta_bin.first, eta_bin.second), "");
            for (size_t idr = 0; idr < dr_bins.size(); ++idr) {
                const auto& dr_bin = dr_bins.at(idr);
                leg.AddEntry(lines.at(idr), Form("#DeltaR #in [%.1f, %.1f]", dr_bin.first, dr_bin.second), "lep");
            }
            leg.Draw();

            for (TH1D* h : lines) {
                delete h;
            }
        }

        std::string full_path = output_dir + "/" + png_name;
        c.SaveAs(full_path.c_str());
        std::cout << "[INFO] Saved: " << full_path << std::endl;
    }
};
