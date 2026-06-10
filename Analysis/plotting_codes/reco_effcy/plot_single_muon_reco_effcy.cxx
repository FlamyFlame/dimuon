#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <ROOT/RDataFrame.hxx>
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <stdexcept>
#include "../../RDFBasedHistFilling/CommonEffcyConfig.h"
#include "../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../MuonObjectsParamsAndHelpers/PbPbBaseClass.h"

struct CtrBinDef {
    int lo, hi;
    std::string suffix;
    std::string label;
};

static std::vector<CtrBinDef> BuildCtrBins(const std::string& ctr_binning_version) {
    struct Helper : public PbPbBaseClass<Helper> {
        friend class PbPbBaseClass<Helper>;
        int RunYear() const { return 24; }
        using PbPbBaseClass<Helper>::InitializePbPb;
        using PbPbBaseClass<Helper>::nCtrBins;
        using PbPbBaseClass<Helper>::ctr_bin_edges;
        using PbPbBaseClass<Helper>::ctr_bins;
    };
    Helper h;
    h.ctr_binning_version = ctr_binning_version;
    h.InitializePbPb();

    std::vector<CtrBinDef> bins;
    bins.push_back({0, 80, "_ctr_inclusive", "0-80%"});
    for (int i = 0; i < h.nCtrBins; ++i)
        bins.push_back({h.ctr_bin_edges[i], h.ctr_bin_edges[i+1],
                        h.ctr_bins[i],
                        std::to_string(h.ctr_bin_edges[i]) + "-" + std::to_string(h.ctr_bin_edges[i+1]) + "%"});
    return bins;
}

static const int colors[] = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1,
                              kOrange+1, kCyan+2, kViolet+1, kTeal+2};

void plot_single_muon_reco_effcy(
    const std::string& mode = "pp",
    const std::string& ctr_binning_version = "default")
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const bool is_overlay = (mode == "overlay");

    std::string data_dir, plot_dir_base, input_file, sample_label;
    if (is_overlay) {
        data_dir   = FullSimSampleInputDir(FullSimSampleType::hijing);
        input_file = data_dir + "muon_pairs_pythia_fullsim_hijing_overlay_pp24_no_data_resonance_cuts_single_muon.root";
        plot_dir_base = data_dir + "plots/hijing_overlay_pp24_single_muon_reco_effcy/";
        sample_label = "Pythia fullsim HIJING overlay";
    } else {
        data_dir   = FullSimSampleInputDir(FullSimSampleType::pp);
        input_file = data_dir + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_single_muon.root";
        plot_dir_base = data_dir + "plots/pp24_single_muon_reco_effcy/";
        sample_label = "Pythia fullsim pp24";
    }

    TFile* f = TFile::Open(input_file.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open " << input_file << "\n";
        return;
    }
    TTree* tree = dynamic_cast<TTree*>(f->Get("muon_tree"));
    if (!tree) {
        std::cerr << "ERROR: muon_tree not found in " << input_file << "\n";
        f->Close(); return;
    }
    std::cout << "Input: " << input_file << " (" << tree->GetEntries() << " entries)\n";
    f->Close();
    delete f;

    static const CommonEffcyConfig ecfg;
    const QEtaBinning& qeta_bins = ecfg.q_eta_proj_ranges_coarse_incl_gap;
    const int nQEta = (int)qeta_bins.size();

    const std::vector<double> pt_edges = {4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 50, 80, 120};
    const int nPtBins = (int)pt_edges.size() - 1;

    std::vector<CtrBinDef> ctr_bins;
    if (is_overlay) {
        ctr_bins = BuildCtrBins(ctr_binning_version);
    } else {
        ctr_bins.push_back({0, 100, "", "pp24"});
    }

    ROOT::RDataFrame df("muon_tree", input_file);

    auto df_base = df.Define("q_eta", "truth_charge * truth_eta");

    ROOT::RDF::RNode df_with_ctr = df_base;

    gSystem->mkdir(plot_dir_base.c_str(), kTRUE);

    for (const auto& ctr : ctr_bins) {
        ROOT::RDF::RNode df_ctr = df_with_ctr;
        if (is_overlay) {
            int lo = ctr.lo, hi = ctr.hi;
            df_ctr = df_with_ctr.Filter([lo, hi](int c){ return c >= lo && c < hi; }, {"ev_centrality"});
        }

        std::vector<TH1D*> h_eff_per_qeta;
        std::vector<std::string> qeta_labels;

        for (int iq = 0; iq < nQEta; ++iq) {
            float qeta_lo = qeta_bins[iq].first;
            float qeta_hi = qeta_bins[iq].second;

            auto df_qeta = df_ctr.Filter(
                [qeta_lo, qeta_hi](float qe){ return qe >= qeta_lo && qe < qeta_hi; },
                {"q_eta"});

            std::string tag = ctr.suffix + "_qeta" + std::to_string(iq);

            auto h_denom = df_qeta.Histo1D(
                {"h_denom", "", nPtBins, pt_edges.data()}, "truth_pt", "ev_weight");
            auto h_num = df_qeta.Filter("pass_medium").Histo1D(
                {"h_num", "", nPtBins, pt_edges.data()}, "truth_pt", "ev_weight");

            TH1D* hd = (TH1D*)h_denom->Clone(("h_denom_" + tag).c_str());
            TH1D* hn = (TH1D*)h_num->Clone(("h_num_" + tag).c_str());
            hd->SetDirectory(nullptr);
            hn->SetDirectory(nullptr);

            TH1D* h_eff = (TH1D*)hn->Clone(("h_eff_" + tag).c_str());
            h_eff->SetDirectory(nullptr);
            h_eff->Divide(hn, hd, 1.0, 1.0, "B");

            for (int ib = 1; ib <= h_eff->GetNbinsX(); ++ib) {
                double den_val = hd->GetBinContent(ib);
                double den_err = hd->GetBinError(ib);
                if (den_err <= 0 || (den_val / den_err) * (den_val / den_err) < 10.0) {
                    h_eff->SetBinContent(ib, 0);
                    h_eff->SetBinError(ib, 0);
                }
            }

            h_eff_per_qeta.push_back(h_eff);
            char buf[64];
            snprintf(buf, sizeof(buf), "%.1f < q#eta < %.1f", qeta_lo, qeta_hi);
            qeta_labels.push_back(buf);

            delete hd; delete hn;
        }

        TCanvas c("c", "c", 1100, 800);
        c.SetTicks(1,1);
        c.SetLogx(true);

        double ymin = 0.0, ymax = 1.15;
        TH1D* h_frame = new TH1D("frame", "", nPtBins, pt_edges.data());
        h_frame->SetMinimum(ymin);
        h_frame->SetMaximum(ymax);
        h_frame->GetXaxis()->SetTitle("Truth p_{T} [GeV]");
        h_frame->GetYaxis()->SetTitle("Single-muon reco efficiency (medium WP)");
        h_frame->GetXaxis()->SetTitleSize(0.045);
        h_frame->GetYaxis()->SetTitleSize(0.045);
        h_frame->GetXaxis()->SetLabelSize(0.038);
        h_frame->GetYaxis()->SetLabelSize(0.038);
        h_frame->GetXaxis()->SetMoreLogLabels(true);
        h_frame->GetXaxis()->SetNoExponent(true);
        h_frame->Draw();

        TLegend leg(0.55, 0.15, 0.90, 0.15 + 0.035 * nQEta);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.028);

        for (int iq = 0; iq < nQEta; ++iq) {
            int ci = colors[iq % 9];
            h_eff_per_qeta[iq]->SetLineColor(ci);
            h_eff_per_qeta[iq]->SetMarkerColor(ci);
            h_eff_per_qeta[iq]->SetMarkerStyle(20 + iq % 10);
            h_eff_per_qeta[iq]->SetMarkerSize(0.8);
            h_eff_per_qeta[iq]->SetLineWidth(2);
            h_eff_per_qeta[iq]->Draw("PE1 SAME");
            leg.AddEntry(h_eff_per_qeta[iq], qeta_labels[iq].c_str(), "lpe");
        }
        leg.Draw("SAME");

        TLatex lat;
        lat.SetNDC();
        lat.SetTextSize(0.040);
        lat.DrawLatex(0.16, 0.92, sample_label.c_str());
        if (is_overlay) {
            lat.SetTextSize(0.036);
            lat.DrawLatex(0.16, 0.87, ("Centrality: " + ctr.label).c_str());
        }

        std::string outname = plot_dir_base + "single_muon_reco_effcy_vs_pt" + ctr.suffix + ".png";
        c.SaveAs(outname.c_str());
        std::cout << "  Saved: " << outname << "\n";

        delete h_frame;
        for (auto* h : h_eff_per_qeta) delete h;
    }

    std::cout << "Done.\n";
}
