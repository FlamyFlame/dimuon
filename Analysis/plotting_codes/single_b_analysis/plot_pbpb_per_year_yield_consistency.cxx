// Per-year CONSISTENCY sanity plot for the combined Pb+Pb cross-section.
//
// For each Run-3 Pb+Pb year found on disk, the signal-region OS pair COUNTS per unit sampled
// luminosity, N_OS / L_int [pairs / nb^-1], as a function of pair pT, one panel per
// ParamsSet::ctrbins centrality class, and below it the ratio of each year to the
// luminosity-weighted all-year mean (= the combined result's raw input).  The years share the
// trigger, the detector and the centrality calibration (D6: PbPb2023 FCal thresholds for all),
// so the ratios should be flat and near 1; a centrality-ordered offset is the signature of a
// per-year FCal-scale / centrality-calibration difference, a pT-dependent one of a per-year
// trigger-efficiency difference.  Raw counts on purpose -- no per-year efficiency correction is
// applied, so what is compared is what the detector delivered.
//
// This is a SANITY plot (per-year consistency), not a per-year cross-section output family: the
// Pb+Pb crossx is always the combined result (feedback_pbpb_crossx_combined).
//
// Inputs : <year>/histograms_real_pairs_pbpb_20YY_<trig>_no_trg_plots_nominal.root
//          h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_<ctr>_counts (X = pair pT, Y = pair eta)
// Binning: the histogram's own pair-pT axis (ParamsSet::pT_bins_150 as booked by the producer),
//          projected as-is -- never retyped here.  Centrality classes from ParamsSet::ctrbins.
// Lumi   : Utilities/PbPbSampledLumi.h (the R_AA luminosity, the same the crossx uses).
// Output : plots/single_b_analysis/pbpb_<years>_combined/sanity/per_year_yield_consistency.png
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"

#include "../../MuonObjectsParamsAndHelpers/DatasetTriggerMap.h"
#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../Utilities/PbPbSampledLumi.h"

void plot_pbpb_per_year_yield_consistency()
{
    const std::string data_base  = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
    const std::string plots_base = data_base + "plots/single_b_analysis/";

    // ---- discover the years (probe-skip-relabel, as plot_single_b_crossx_pbpb does) ----
    struct Year { int yr; double lumi; TFile* f; int color; int marker; };
    std::vector<Year> years;
    const int colors[]  = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1};
    const int markers[] = {20, 21, 22, 23};
    std::string year_tag;
    for (int yr : {23, 24, 25, 26}) {
        const std::string& trig = DatasetTriggerMap::GetTrigger(yr, "PbPb");
        const std::string path = data_base + "pbpb_20" + std::to_string(yr)
            + "/histograms_real_pairs_pbpb_20" + std::to_string(yr) + "_" + trig
            + "_no_trg_plots_nominal.root";
        if (gSystem->AccessPathName(path.c_str())) {
            std::cout << "[INFO] PbPb 20" << yr << ": no crossx file -- skipped." << std::endl;
            continue;
        }
        TFile* f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + path);
        const size_t i = years.size();
        years.push_back({yr, PbPbMu4SampledLumiNb(yr), f, colors[i % 4], markers[i % 4]});
        year_tag += (year_tag.empty() ? "" : "_") + std::to_string(yr);
    }
    if (years.size() < 2)
        throw std::runtime_error("per-year consistency needs at least two Pb+Pb years on disk");

    const std::string out_dir = plots_base + "pbpb_" + year_tag + "_combined/sanity";
    gSystem->mkdir(out_dir.c_str(), true);

    // ---- centrality classes from ParamsSet::ctrbins ----
    const auto& cb = ParamsSet::ctrbins;
    const int n_ctr = static_cast<int>(cb.size()) - 1;

    gStyle->SetOptStat(0);
    // nrows >= ncols, nrows ~ sqrt(N): 6 classes x 2 (yield + ratio) -> 3 columns x 4 rows.
    const int ncols = 3, nrows = 4;
    TCanvas c("c_per_year_yield", "", 600 * ncols, 500 * nrows);
    c.Divide(ncols, nrows);

    std::vector<TH1D*> keep;  // owned until SaveAs
    TLegend* leg = nullptr;
    for (int ic = 0; ic < n_ctr; ++ic) {
        const std::string ctr = "ctr" + std::to_string(cb[ic]) + "_" + std::to_string(cb[ic + 1]);
        const std::string hname =
            "h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr + "_counts";

        // per-year N_OS / L  vs pair pT, on the histogram's own pair-pT axis
        std::vector<TH1D*> per_year;
        double lumi_sum = 0.;
        TH1D* sum = nullptr;  // lumi-weighted mean numerator: sum_y N_y
        for (auto& y : years) {
            TH2D* h2 = dynamic_cast<TH2D*>(y.f->Get(hname.c_str()));
            if (!h2) throw std::runtime_error("PbPb 20" + std::to_string(y.yr) + ": missing " + hname);
            TH1D* h = h2->ProjectionX(Form("h_yield_%s_%d", ctr.c_str(), y.yr));
            h->SetDirectory(nullptr);
            if (!sum) { sum = static_cast<TH1D*>(h->Clone(Form("h_sum_%s", ctr.c_str()))); sum->SetDirectory(nullptr); }
            else sum->Add(h);
            lumi_sum += y.lumi;
            h->Scale(1.0 / y.lumi);
            h->SetLineColor(y.color); h->SetMarkerColor(y.color); h->SetMarkerStyle(y.marker);
            h->SetMarkerSize(1.0);
            per_year.push_back(h);
            keep.push_back(h);
        }
        // combined = (sum_y N_y) / (sum_y L_y): the raw input of the combined cross-section
        sum->Scale(1.0 / lumi_sum);
        keep.push_back(sum);

        // yield panel
        c.cd(ic + 1);
        gPad->SetLogx(); gPad->SetLogy();
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        double ymax = 0.;
        for (auto* h : per_year) ymax = std::max(ymax, h->GetBinContent(h->GetMaximumBin()));
        per_year[0]->SetTitle(Form("Pb+Pb %d-%d%%;p_{T}^{#mu#mu} [GeV];N_{OS}^{signal region} / L_{int}  [pairs / nb^{-1}]",
                                   cb[ic], cb[ic + 1]));
        per_year[0]->SetMaximum(ymax * 3.);
        per_year[0]->SetMinimum(0.05);
        per_year[0]->Draw("E1");
        for (size_t i = 1; i < per_year.size(); ++i) per_year[i]->Draw("E1 same");
        if (ic == 0) {
            leg = new TLegend(0.55, 0.62, 0.88, 0.88);
            leg->SetBorderSize(0); leg->SetFillStyle(0);
            for (size_t i = 0; i < years.size(); ++i)
                leg->AddEntry(per_year[i], Form("Pb+Pb 20%d  (%.3f nb^{-1})", years[i].yr, years[i].lumi), "lp");
            leg->Draw();
        }

        // ratio panel: year / lumi-weighted all-year mean
        c.cd(n_ctr + ic + 1);
        gPad->SetLogx();
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
        bool first = true;
        for (auto* h : per_year) {
            TH1D* r = static_cast<TH1D*>(h->Clone(Form("%s_ratio", h->GetName())));
            r->SetDirectory(nullptr);
            r->Divide(sum);
            r->SetTitle(Form("Pb+Pb %d-%d%%;p_{T}^{#mu#mu} [GeV];year / lumi-weighted mean", cb[ic], cb[ic + 1]));
            r->SetMinimum(0.5); r->SetMaximum(1.5);
            r->Draw(first ? "E1" : "E1 same");
            first = false;
            keep.push_back(r);
        }
        TLine* l = new TLine(per_year[0]->GetXaxis()->GetXmin(), 1., per_year[0]->GetXaxis()->GetXmax(), 1.);
        l->SetLineStyle(2); l->Draw();
    }

    const std::string out = out_dir + "/per_year_yield_consistency.png";
    c.SaveAs(out.c_str());
    std::cout << "[INFO] Saved: " << out << std::endl;

    // integrated numbers, for the log
    std::cout << "\nSignal-region OS pairs per nb^-1 (integrated over pair pT):" << std::endl;
    for (int ic = 0; ic < n_ctr; ++ic) {
        const std::string ctr = "ctr" + std::to_string(cb[ic]) + "_" + std::to_string(cb[ic + 1]);
        std::cout << Form("  %-9s", ctr.c_str());
        for (auto& y : years) {
            TH2D* h2 = dynamic_cast<TH2D*>(y.f->Get(("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_" + ctr + "_counts").c_str()));
            std::cout << Form("  20%d: %7.0f", y.yr, h2->Integral() / y.lumi);
        }
        std::cout << std::endl;
    }
    for (auto& y : years) y.f->Close();
}
