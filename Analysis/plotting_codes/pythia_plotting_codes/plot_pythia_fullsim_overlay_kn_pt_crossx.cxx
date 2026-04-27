// Differential cross-section vs pair pT per pTHat (kn) range — Pythia fullsim HIJING overlay.
// Mirrors plot_pythia_fullsim_kn_pt_crossx.cxx; scale_factors = N_full/N_test(10k):
//   kn0-5: {160, 320, 200, 60, 40, 40}
// Selection: muon_pair_tree_kin*_sign2 with from_same_b.
// Two plots: truth_pair_pt and reco pair_pt (pair_pass_medium additionally required).

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <THStack.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TPad.h>
#include <cmath>
#include <string>
#include <vector>
#include <array>

static const std::string kInputFile =
    "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/"
    "muon_pairs_pythia_fullsim_hijing_overlay_pp24_no_data_resonance_cuts.root";
static const std::string kOutputDir =
    "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/plots/";

static const int    kNkn = 6;
static const std::array<std::string, 6> kKnLabels = {
    "#hat{p}_{T} 8-14 GeV",
    "#hat{p}_{T} 14-24 GeV",
    "#hat{p}_{T} 24-40 GeV",
    "#hat{p}_{T} 40-70 GeV",
    "#hat{p}_{T} 70-125 GeV",
    "#hat{p}_{T} 125-300 GeV"
};
static const std::array<std::string, 6> kKnShort = {
    "8-14", "14-24", "24-40", "40-70", "70-125", "125-300"
};
static const std::array<int, 6> kColors = {
    kRed+1, kOrange+1, kGreen+2, kCyan+2, kBlue+1, kViolet+1
};
// N_full / N_test(10k): {1600k,3200k,2000k,600k,400k,400k} / 10k
static const std::array<double, 6> kScaleFactors = {160., 320., 200., 60., 40., 40.};

static const std::array<std::string, 2> kVarNames   = { "truth_pair_pt", "pair_pt" };
static const std::array<std::string, 2> kVarTitles  = {
    "truth p_{T}^{pair} [GeV]", "reco p_{T}^{pair} [GeV]"
};
static const std::array<std::string, 2> kFilters = {
    "from_same_b",
    "from_same_b && pair_pass_medium"
};
static const std::array<std::string, 2> kLabels = {
    "truth p_{T}^{pair}", "reco p_{T}^{pair}"
};

static std::vector<double> LogEdges(int nbins, double xmin, double xmax) {
    std::vector<double> e(nbins + 1);
    const double lmin = std::log(xmin), lmax = std::log(xmax);
    for (int i = 0; i <= nbins; i++)
        e[i] = std::exp(lmin + i * (lmax - lmin) / nbins);
    return e;
}

// ---------------------------------------------------------------------------
void overlay_plot_impl(int nbins_arg, double xmax_arg, const std::string& suffix) {

    gSystem->mkdir(kOutputDir.c_str(), true);
    const auto edges = LogEdges(nbins_arg, 8., xmax_arg);

    const std::array<std::string, 2> out_names = {
        "truth_pair_pt_kn" + suffix, "reco_pair_pt_kn" + suffix
    };

    for (int ivar = 0; ivar < 2; ivar++) {
        std::vector<TH1D*> hists(kNkn);
        for (int ikn = 0; ikn < kNkn; ikn++) {
            const std::string tree = "muon_pair_tree_kin" + std::to_string(ikn) + "_sign2";
            ROOT::RDataFrame df(tree, kInputFile);
            auto hptr = df.Filter(kFilters[ivar])
                          .Histo1D(ROOT::RDF::TH1DModel{
                              ("h_" + kVarNames[ivar] + "_kn" + std::to_string(ikn)).c_str(), "",
                              nbins_arg, edges.data()
                          }, kVarNames[ivar], "weight");
            hists[ikn] = (TH1D*)hptr->Clone();
            hists[ikn]->SetDirectory(nullptr);
            hists[ikn]->Scale(1., "width");
            hists[ikn]->SetLineColor(kColors[ikn]);
            hists[ikn]->SetMarkerColor(kColors[ikn]);
            hists[ikn]->SetMarkerStyle(20);
            hists[ikn]->SetMarkerSize(0.7);
            hists[ikn]->SetLineWidth(1);
            hists[ikn]->SetFillColor(kColors[ikn]);
            hists[ikn]->SetFillStyle(1001);
        }

        double ymax = 0., ymin = 1e30;
        for (auto* h : hists) {
            ymax = std::max(ymax, h->GetMaximum());
            for (int ib = 1; ib <= h->GetNbinsX(); ib++)
                if (h->GetBinContent(ib) > 0) ymin = std::min(ymin, h->GetBinContent(ib));
        }
        if (ymin > 1e29) ymin = 1e-12;

        TCanvas* c = new TCanvas(out_names[ivar].c_str(), "", 1400, 600);
        c->Divide(2, 1);

        c->cd(1);
        gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14); gPad->SetLogx(); gPad->SetLogy();
        for (int ikn = 0; ikn < kNkn; ikn++) {
            auto* h = hists[ikn];
            h->GetXaxis()->SetTitle(kVarTitles[ivar].c_str());
            h->GetYaxis()->SetTitle("d#sigma/dp_{T} [#mub/GeV]");
            h->GetXaxis()->SetRangeUser(8., xmax_arg);
            h->GetYaxis()->SetRangeUser(ymin * 0.3, ymax * 5.);
            h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
            h->GetXaxis()->SetTitleOffset(1.1); h->GetYaxis()->SetTitleOffset(1.5);
            if (ikn == 0) h->Draw("E"); else h->Draw("E same");
        }
        TLegend* leg1 = new TLegend(0.67, 0.53, 1.08, 0.92);
        leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextSize(0.034);
        for (int ikn = 0; ikn < kNkn; ikn++)
            leg1->AddEntry(hists[ikn], kKnLabels[ikn].c_str(), "lep");
        leg1->Draw();
        TLatex lat1; lat1.SetNDC(); lat1.SetTextSize(0.038);
        lat1.DrawLatex(0.17, 0.92, ("Pythia fullsim hijing overlay pp24, single-b, " + kLabels[ivar]).c_str());

        c->cd(2);
        gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14); gPad->SetLogx(); gPad->SetLogy();
        THStack* hs = new THStack(("hs_" + out_names[ivar]).c_str(), "");
        for (int ikn = 0; ikn < kNkn; ikn++) hs->Add(hists[ikn]);
        hs->Draw("hist");
        hs->GetXaxis()->SetTitle(kVarTitles[ivar].c_str());
        hs->GetYaxis()->SetTitle("d#sigma/dp_{T} [#mub/GeV]");
        hs->GetXaxis()->SetRangeUser(8., xmax_arg);
        hs->GetXaxis()->SetTitleSize(0.05); hs->GetYaxis()->SetTitleSize(0.05);
        hs->GetXaxis()->SetTitleOffset(1.1); hs->GetYaxis()->SetTitleOffset(1.5);
        TLegend* leg2 = new TLegend(0.67, 0.53, 1.08, 0.92);
        leg2->SetBorderSize(0); leg2->SetFillStyle(0); leg2->SetTextSize(0.034);
        for (int ikn = 0; ikn < kNkn; ikn++)
            leg2->AddEntry(hists[ikn], kKnLabels[ikn].c_str(), "f");
        leg2->Draw();
        TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(0.038);
        lat2.DrawLatex(0.17, 0.92, ("Pythia fullsim hijing overlay pp24, single-b, " + kLabels[ivar]).c_str());

        c->SaveAs((kOutputDir + out_names[ivar] + ".png").c_str());
        for (auto* h : hists) delete h;
        delete hs; delete c;
    }
}

// ---------------------------------------------------------------------------
void overlay_plot_stat_error_forecast(int nbins_arg, double xmax_arg, const std::string& suffix,
                                       const std::string& subdir = "") {

    gSystem->mkdir((kOutputDir + subdir).c_str(), true);
    const auto edges = LogEdges(nbins_arg, 8., xmax_arg);

    const std::array<std::string, 2> out_names = {
        "truth_pair_pt_kn_stat_error" + suffix,
        "reco_pair_pt_kn_stat_error"  + suffix
    };

    for (int ivar = 0; ivar < 2; ivar++) {
        std::vector<TH1D*> hists(kNkn);
        for (int ikn = 0; ikn < kNkn; ikn++) {
            const std::string tree = "muon_pair_tree_kin" + std::to_string(ikn) + "_sign2";
            ROOT::RDataFrame df(tree, kInputFile);
            auto hptr = df.Filter(kFilters[ivar])
                          .Histo1D(ROOT::RDF::TH1DModel{
                              ("hse_" + kVarNames[ivar] + "_kn" + std::to_string(ikn)).c_str(), "",
                              nbins_arg, edges.data()
                          }, kVarNames[ivar], "weight");
            hists[ikn] = (TH1D*)hptr->Clone();
            hists[ikn]->SetDirectory(nullptr);
            hists[ikn]->Scale(1., "width");
            for (int ib = 1; ib <= hists[ikn]->GetNbinsX(); ib++)
                hists[ikn]->SetBinError(ib, hists[ikn]->GetBinError(ib) / std::sqrt(kScaleFactors[ikn]));
            hists[ikn]->SetLineColor(kColors[ikn]); hists[ikn]->SetMarkerColor(kColors[ikn]);
            hists[ikn]->SetMarkerStyle(20); hists[ikn]->SetMarkerSize(0.7); hists[ikn]->SetLineWidth(1);
        }

        std::vector<TH1D*> herr(kNkn);
        for (int ikn = 0; ikn < kNkn; ikn++) {
            herr[ikn] = (TH1D*)hists[ikn]->Clone(
                ("herr_" + kVarNames[ivar] + "_kn" + std::to_string(ikn) + suffix).c_str());
            herr[ikn]->SetDirectory(nullptr); herr[ikn]->Reset();
            for (int ib = 1; ib <= hists[ikn]->GetNbinsX(); ib++) {
                const double cv = hists[ikn]->GetBinContent(ib);
                const double e  = hists[ikn]->GetBinError(ib);
                herr[ikn]->SetBinContent(ib, cv > 0 ? e / cv : 0.);
                herr[ikn]->SetBinError(ib, 0.);
            }
            herr[ikn]->SetLineColor(kColors[ikn]); herr[ikn]->SetLineWidth(2);
        }

        double ymax_d = 0., ymin_d = 1e30, ymax_e = 0., ymin_e = 1e30;
        for (auto* h : hists) {
            ymax_d = std::max(ymax_d, h->GetMaximum());
            for (int ib = 1; ib <= h->GetNbinsX(); ib++)
                if (h->GetBinContent(ib) > 0) ymin_d = std::min(ymin_d, h->GetBinContent(ib));
        }
        for (auto* h : herr) {
            ymax_e = std::max(ymax_e, h->GetMaximum());
            for (int ib = 1; ib <= h->GetNbinsX(); ib++)
                if (h->GetBinContent(ib) > 0) ymin_e = std::min(ymin_e, h->GetBinContent(ib));
        }
        if (ymin_d > 1e29) ymin_d = 1e-12;
        if (ymin_e > 1e29) ymin_e = 1e-12;

        TCanvas* c = new TCanvas(out_names[ivar].c_str(), "", 1400, 600);
        c->Divide(2, 1);

        c->cd(1);
        gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14); gPad->SetLogx(); gPad->SetLogy();
        for (int ikn = 0; ikn < kNkn; ikn++) {
            auto* h = hists[ikn];
            h->GetXaxis()->SetTitle(kVarTitles[ivar].c_str());
            h->GetYaxis()->SetTitle("d#sigma/dp_{T} [#mub/GeV]");
            h->GetXaxis()->SetRangeUser(8., xmax_arg);
            h->GetYaxis()->SetRangeUser(ymin_d * 0.3, ymax_d * 5.);
            h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
            h->GetXaxis()->SetTitleOffset(1.1); h->GetYaxis()->SetTitleOffset(1.5);
            if (ikn == 0) h->Draw("E"); else h->Draw("E same");
        }
        TLegend* leg1 = new TLegend(0.67, 0.53, 1.08, 0.92);
        leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextSize(0.034);
        for (int ikn = 0; ikn < kNkn; ikn++)
            leg1->AddEntry(hists[ikn], kKnLabels[ikn].c_str(), "lep");
        leg1->Draw();
        TLatex lat1; lat1.SetNDC(); lat1.SetTextSize(0.038);
        lat1.DrawLatex(0.17, 0.92, ("Pythia fullsim hijing overlay pp24, single-b, " + kLabels[ivar]).c_str());

        c->cd(2);
        gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14); gPad->SetLogx(); gPad->SetLogy();
        for (int ikn = 0; ikn < kNkn; ikn++) {
            auto* h = herr[ikn];
            h->GetXaxis()->SetTitle(kVarTitles[ivar].c_str());
            h->GetYaxis()->SetTitle("Full-sample rel. stat. error 1/#sqrt{N_{full}}");
            h->GetXaxis()->SetRangeUser(8., xmax_arg);
            h->GetYaxis()->SetRangeUser(ymin_e * 0.3, ymax_e * 5.);
            h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
            h->GetXaxis()->SetTitleOffset(1.1); h->GetYaxis()->SetTitleOffset(1.5);
            if (ikn == 0) h->Draw("hist"); else h->Draw("hist same");
        }
        TLegend* leg2 = new TLegend(0.67, 0.53, 1.08, 0.92);
        leg2->SetBorderSize(0); leg2->SetFillStyle(0); leg2->SetTextSize(0.034);
        for (int ikn = 0; ikn < kNkn; ikn++)
            leg2->AddEntry(herr[ikn], kKnLabels[ikn].c_str(), "l");
        leg2->Draw();
        TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(0.038);
        lat2.DrawLatex(0.17, 0.92, ("Rel. stat. error (full sample), " + kLabels[ivar]).c_str());

        c->SaveAs((kOutputDir + subdir + out_names[ivar] + ".png").c_str());
        for (auto* h : hists) delete h;
        for (auto* h : herr)  delete h;
        delete c;
    }
}

// ---------------------------------------------------------------------------
void overlay_plot_err_fraction_map(int nbins_arg, double xmax_arg, const std::string& suffix,
                                    const std::string& subdir = "") {

    gSystem->mkdir((kOutputDir + subdir).c_str(), true);
    const auto edges = LogEdges(nbins_arg, 8., xmax_arg);

    const std::array<std::string, 2> out_names = {
        "truth_pair_pt_kn_err_frac" + suffix,
        "reco_pair_pt_kn_err_frac"  + suffix
    };

    for (int ivar = 0; ivar < 2; ivar++) {
        std::vector<TH1D*> hists(kNkn);
        for (int ikn = 0; ikn < kNkn; ikn++) {
            const std::string tree = "muon_pair_tree_kin" + std::to_string(ikn) + "_sign2";
            ROOT::RDataFrame df(tree, kInputFile);
            auto hptr = df.Filter(kFilters[ivar])
                          .Histo1D(ROOT::RDF::TH1DModel{
                              ("hef_" + kVarNames[ivar] + "_kn" + std::to_string(ikn)).c_str(), "",
                              nbins_arg, edges.data()
                          }, kVarNames[ivar], "weight");
            hists[ikn] = (TH1D*)hptr->Clone();
            hists[ikn]->SetDirectory(nullptr);
            hists[ikn]->Scale(1., "width");
            for (int ib = 1; ib <= hists[ikn]->GetNbinsX(); ib++)
                hists[ikn]->SetBinError(ib, hists[ikn]->GetBinError(ib) / std::sqrt(kScaleFactors[ikn]));
            hists[ikn]->SetLineColor(kColors[ikn]); hists[ikn]->SetMarkerColor(kColors[ikn]);
            hists[ikn]->SetMarkerStyle(20); hists[ikn]->SetMarkerSize(0.7); hists[ikn]->SetLineWidth(1);
        }

        TH2D* h2 = new TH2D(("h2ef_" + out_names[ivar]).c_str(), "",
                             nbins_arg, edges.data(), kNkn, -0.5, kNkn - 0.5);
        h2->SetDirectory(nullptr);
        for (int ib = 1; ib <= nbins_arg; ib++) {
            double tot2 = 0.;
            for (int ikn = 0; ikn < kNkn; ikn++) { double e = hists[ikn]->GetBinError(ib); tot2 += e * e; }
            const double tot = std::sqrt(tot2);
            for (int ikn = 0; ikn < kNkn; ikn++) {
                double e = hists[ikn]->GetBinError(ib);
                h2->SetBinContent(ib, ikn + 1, tot > 0 ? e / tot : 0.);
            }
        }
        for (int ikn = 0; ikn < kNkn; ikn++)
            h2->GetYaxis()->SetBinLabel(ikn + 1, kKnShort[ikn].c_str());

        double ymax_d = 0., ymin_d = 1e30;
        for (auto* h : hists) {
            ymax_d = std::max(ymax_d, h->GetMaximum());
            for (int ib = 1; ib <= h->GetNbinsX(); ib++)
                if (h->GetBinContent(ib) > 0) ymin_d = std::min(ymin_d, h->GetBinContent(ib));
        }
        if (ymin_d > 1e29) ymin_d = 1e-12;

        TCanvas* c = new TCanvas(out_names[ivar].c_str(), "", 1400, 600);
        c->Divide(2, 1);

        c->cd(1);
        gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14); gPad->SetLogx(); gPad->SetLogy();
        for (int ikn = 0; ikn < kNkn; ikn++) {
            auto* h = hists[ikn];
            h->GetXaxis()->SetTitle(kVarTitles[ivar].c_str());
            h->GetYaxis()->SetTitle("d#sigma/dp_{T} [#mub/GeV]");
            h->GetXaxis()->SetRangeUser(8., xmax_arg);
            h->GetYaxis()->SetRangeUser(ymin_d * 0.3, ymax_d * 5.);
            h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
            h->GetXaxis()->SetTitleOffset(1.1); h->GetYaxis()->SetTitleOffset(1.5);
            if (ikn == 0) h->Draw("E"); else h->Draw("E same");
        }
        TLegend* leg1 = new TLegend(0.67, 0.53, 1.08, 0.92);
        leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextSize(0.034);
        for (int ikn = 0; ikn < kNkn; ikn++)
            leg1->AddEntry(hists[ikn], kKnLabels[ikn].c_str(), "lep");
        leg1->Draw();
        TLatex lat1; lat1.SetNDC(); lat1.SetTextSize(0.038);
        lat1.DrawLatex(0.17, 0.92, ("Pythia fullsim hijing overlay pp24, single-b, " + kLabels[ivar]).c_str());

        c->cd(2);
        gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.16);
        gPad->SetBottomMargin(0.14); gPad->SetLogx();
        gStyle->SetPalette(kViridis);
        h2->GetXaxis()->SetTitle(kVarTitles[ivar].c_str());
        h2->GetYaxis()->SetTitle("#hat{p}_{T} range [GeV]");
        h2->GetZaxis()->SetTitle("#sigma_{err,kn} / #sigma_{err,total}");
        h2->GetXaxis()->SetRangeUser(8., xmax_arg);
        h2->SetMinimum(0.); h2->SetMaximum(1.);
        h2->GetXaxis()->SetTitleSize(0.05); h2->GetYaxis()->SetTitleSize(0.05);
        h2->GetZaxis()->SetTitleSize(0.045);
        h2->GetXaxis()->SetTitleOffset(1.1); h2->GetYaxis()->SetTitleOffset(1.2);
        h2->GetZaxis()->SetTitleOffset(1.3); h2->GetYaxis()->SetLabelSize(0.05);
        h2->Draw("colz");
        TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(0.038);
        lat2.DrawLatex(0.15, 0.92, ("Full-sample error fraction, " + kLabels[ivar]).c_str());

        c->SaveAs((kOutputDir + subdir + out_names[ivar] + ".png").c_str());
        for (auto* h : hists) delete h;
        delete h2; delete c;
    }
}

// ---------------------------------------------------------------------------
void overlay_plot_err_ratio_map(int nbins_arg, double xmax_arg, const std::string& suffix,
                                 const std::string& subdir = "") {

    gSystem->mkdir((kOutputDir + subdir).c_str(), true);
    const auto edges = LogEdges(nbins_arg, 8., xmax_arg);

    const std::array<std::string, 2> out_names = {
        "truth_pair_pt_kn_err_ratio" + suffix,
        "reco_pair_pt_kn_err_ratio"  + suffix
    };

    for (int ivar = 0; ivar < 2; ivar++) {
        std::vector<TH1D*> hists(kNkn);
        for (int ikn = 0; ikn < kNkn; ikn++) {
            const std::string tree = "muon_pair_tree_kin" + std::to_string(ikn) + "_sign2";
            ROOT::RDataFrame df(tree, kInputFile);
            auto hptr = df.Filter(kFilters[ivar])
                          .Histo1D(ROOT::RDF::TH1DModel{
                              ("her_" + kVarNames[ivar] + "_kn" + std::to_string(ikn)).c_str(), "",
                              nbins_arg, edges.data()
                          }, kVarNames[ivar], "weight");
            hists[ikn] = (TH1D*)hptr->Clone();
            hists[ikn]->SetDirectory(nullptr);
            hists[ikn]->Scale(1., "width");
            for (int ib = 1; ib <= hists[ikn]->GetNbinsX(); ib++)
                hists[ikn]->SetBinError(ib, hists[ikn]->GetBinError(ib) / std::sqrt(kScaleFactors[ikn]));
            hists[ikn]->SetLineColor(kColors[ikn]); hists[ikn]->SetMarkerColor(kColors[ikn]);
            hists[ikn]->SetMarkerStyle(20); hists[ikn]->SetMarkerSize(0.7); hists[ikn]->SetLineWidth(1);
        }

        TH2D* h2 = new TH2D(("h2er_" + out_names[ivar]).c_str(), "",
                             nbins_arg, edges.data(), kNkn, -0.5, kNkn - 0.5);
        h2->SetDirectory(nullptr);
        for (int ib = 1; ib <= nbins_arg; ib++) {
            double tot_err2 = 0., tot_xsec = 0.;
            for (int ikn = 0; ikn < kNkn; ikn++) {
                double e = hists[ikn]->GetBinError(ib);
                tot_err2 += e * e;
                tot_xsec += hists[ikn]->GetBinContent(ib);
            }
            const double tot_err = std::sqrt(tot_err2);
            for (int ikn = 0; ikn < kNkn; ikn++) {
                double e     = hists[ikn]->GetBinError(ib);
                double cv    = hists[ikn]->GetBinContent(ib);
                double ferr  = (tot_err  > 0) ? e  / tot_err  : 0.;
                double fxsec = (tot_xsec > 0) ? cv / tot_xsec : 0.;
                h2->SetBinContent(ib, ikn + 1, fxsec > 0 ? ferr / fxsec : 0.);
            }
        }
        for (int ikn = 0; ikn < kNkn; ikn++)
            h2->GetYaxis()->SetBinLabel(ikn + 1, kKnShort[ikn].c_str());

        double ymax_d = 0., ymin_d = 1e30;
        for (auto* h : hists) {
            ymax_d = std::max(ymax_d, h->GetMaximum());
            for (int ib = 1; ib <= h->GetNbinsX(); ib++)
                if (h->GetBinContent(ib) > 0) ymin_d = std::min(ymin_d, h->GetBinContent(ib));
        }
        if (ymin_d > 1e29) ymin_d = 1e-12;

        TCanvas* c = new TCanvas(out_names[ivar].c_str(), "", 1400, 600);
        c->Divide(2, 1);

        c->cd(1);
        gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14); gPad->SetLogx(); gPad->SetLogy();
        for (int ikn = 0; ikn < kNkn; ikn++) {
            auto* h = hists[ikn];
            h->GetXaxis()->SetTitle(kVarTitles[ivar].c_str());
            h->GetYaxis()->SetTitle("d#sigma/dp_{T} [#mub/GeV]");
            h->GetXaxis()->SetRangeUser(8., xmax_arg);
            h->GetYaxis()->SetRangeUser(ymin_d * 0.3, ymax_d * 5.);
            h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
            h->GetXaxis()->SetTitleOffset(1.1); h->GetYaxis()->SetTitleOffset(1.5);
            if (ikn == 0) h->Draw("E"); else h->Draw("E same");
        }
        TLegend* leg1 = new TLegend(0.67, 0.53, 1.08, 0.92);
        leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextSize(0.034);
        for (int ikn = 0; ikn < kNkn; ikn++)
            leg1->AddEntry(hists[ikn], kKnLabels[ikn].c_str(), "lep");
        leg1->Draw();
        TLatex lat1; lat1.SetNDC(); lat1.SetTextSize(0.038);
        lat1.DrawLatex(0.17, 0.92, ("Pythia fullsim hijing overlay pp24, single-b, " + kLabels[ivar]).c_str());

        c->cd(2);
        gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.16);
        gPad->SetBottomMargin(0.14); gPad->SetLogx();
        gStyle->SetPalette(kViridis);
        h2->GetXaxis()->SetTitle(kVarTitles[ivar].c_str());
        h2->GetYaxis()->SetTitle("#hat{p}_{T} range [GeV]");
        h2->GetZaxis()->SetTitle("(#sigma_{err,kn}/#sigma_{err,tot}) / (#sigma_{kn}/#sigma_{tot})");
        h2->GetXaxis()->SetRangeUser(8., xmax_arg);
        h2->SetMinimum(0.);
        h2->GetXaxis()->SetTitleSize(0.05); h2->GetYaxis()->SetTitleSize(0.05);
        h2->GetZaxis()->SetTitleSize(0.038);
        h2->GetXaxis()->SetTitleOffset(1.1); h2->GetYaxis()->SetTitleOffset(1.2);
        h2->GetZaxis()->SetTitleOffset(1.5); h2->GetYaxis()->SetLabelSize(0.05);
        h2->Draw("colz");
        TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(0.038);
        lat2.DrawLatex(0.15, 0.92, ("Full-sample err/xsec ratio, " + kLabels[ivar]).c_str());

        c->SaveAs((kOutputDir + subdir + out_names[ivar] + ".png").c_str());
        for (auto* h : hists) delete h;
        delete h2; delete c;
    }
}

// ---------------------------------------------------------------------------
void plot_pythia_fullsim_overlay_kn_pt_crossx() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    overlay_plot_impl(20, 120., "");
    overlay_plot_impl(25, 150., "_150GeV");
    overlay_plot_stat_error_forecast(20, 120., "");
    overlay_plot_stat_error_forecast(25, 150., "_150GeV");
    overlay_plot_err_fraction_map(20, 120., "");
    overlay_plot_err_fraction_map(25, 150., "_150GeV");
    overlay_plot_err_ratio_map(20, 120., "");
    overlay_plot_err_ratio_map(25, 150., "_150GeV");
}
