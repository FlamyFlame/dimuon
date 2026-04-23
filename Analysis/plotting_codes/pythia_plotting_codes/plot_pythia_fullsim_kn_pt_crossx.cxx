// Differential cross-section vs pair pT per pTHat (kn) range.
// Selection: muon_pair_tree_kin*_sign2 with from_same_b.
// Two plots: truth_pair_pt and reco pair_pt (pair_pass_medium additionally required).
// Left subplot: markers+errorbars per kn.  Right subplot: stack.
// Entry point calls two binning versions: 20 bins 8-120 GeV and 25 bins 8-150 GeV.

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

void plot_impl(int nbins_arg, double xmax_arg, const std::string& suffix) {

    const std::string input_file =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/"
        "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts.root";
    const std::string output_dir =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/plots/";
    gSystem->mkdir(output_dir.c_str(), true);

    // kn ranges
    const int nkn = 6;
    const std::array<std::string, nkn> kn_labels = {
        "#hat{p}_{T} 8-14 GeV",
        "#hat{p}_{T} 14-24 GeV",
        "#hat{p}_{T} 24-40 GeV",
        "#hat{p}_{T} 40-70 GeV",
        "#hat{p}_{T} 70-125 GeV",
        "#hat{p}_{T} 125-300 GeV"
    };
    const std::array<int, nkn> colors = {
        kRed+1, kOrange+1, kGreen+2, kCyan+2, kBlue+1, kViolet+1
    };

    const int    nbins = nbins_arg;
    const double xmin  = 8., xmax = xmax_arg;
    std::vector<double> edges(nbins + 1);
    const double lmin = std::log(xmin), lmax = std::log(xmax);
    for (int i = 0; i <= nbins; i++)
        edges[i] = std::exp(lmin + i * (lmax - lmin) / nbins);

    // Two variables: [0] truth_pair_pt, [1] pair_pt (reco)
    const std::array<std::string, 2> var_names   = { "truth_pair_pt", "pair_pt" };
    const std::array<std::string, 2> var_titles  = {
        "truth p_{T}^{pair} [GeV]", "reco p_{T}^{pair} [GeV]"
    };
    const std::array<std::string, 2> filter_strs = {
        "from_same_b",
        "from_same_b && pair_pass_medium"
    };
    const std::array<std::string, 2> out_names   = {
        "truth_pair_pt_kn" + suffix, "reco_pair_pt_kn" + suffix
    };

    for (int ivar = 0; ivar < 2; ivar++) {
        const auto& var     = var_names[ivar];
        const auto& xtitle  = var_titles[ivar];
        const auto& filter  = filter_strs[ivar];
        const auto& outname = out_names[ivar];

        std::vector<TH1D*> hists(nkn);
        for (int ikn = 0; ikn < nkn; ikn++) {
            const std::string tree = "muon_pair_tree_kin" + std::to_string(ikn) + "_sign2";
            ROOT::RDataFrame df(tree, input_file);
            auto hptr = df.Filter(filter)
                          .Histo1D(ROOT::RDF::TH1DModel{
                              ("h_" + var + "_kn" + std::to_string(ikn)).c_str(), "",
                              nbins, edges.data()
                          }, var, "weight");
            hists[ikn] = (TH1D*)hptr->Clone();
            hists[ikn]->SetDirectory(nullptr);
            hists[ikn]->Scale(1., "width");   // differential

            hists[ikn]->SetLineColor(colors[ikn]);
            hists[ikn]->SetMarkerColor(colors[ikn]);
            hists[ikn]->SetMarkerStyle(20);
            hists[ikn]->SetMarkerSize(0.7);
            hists[ikn]->SetLineWidth(1);
            hists[ikn]->SetFillColor(colors[ikn]);
            hists[ikn]->SetFillStyle(1001);
        }

        // find global ymax across all kn (for markers panel)
        double ymax = 0.;
        for (auto* h : hists) ymax = std::max(ymax, h->GetMaximum());
        double ymin_nonzero = 1e30;
        for (auto* h : hists)
            for (int ib = 1; ib <= h->GetNbinsX(); ib++)
                if (h->GetBinContent(ib) > 0)
                    ymin_nonzero = std::min(ymin_nonzero, h->GetBinContent(ib));
        if (ymin_nonzero > 1e29) ymin_nonzero = 1e-12;

        TCanvas* c = new TCanvas(outname.c_str(), "", 1400, 600);
        c->Divide(2, 1);

        // --- Left: markers with error bars ---
        c->cd(1);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14);
        gPad->SetLogx();
        gPad->SetLogy();

        for (int ikn = 0; ikn < nkn; ikn++) {
            auto* h = hists[ikn];
            h->GetXaxis()->SetTitle(xtitle.c_str());
            h->GetYaxis()->SetTitle("d#sigma/dp_{T} [#mub/GeV]");
            h->GetXaxis()->SetRangeUser(xmin, xmax);
            h->GetYaxis()->SetRangeUser(ymin_nonzero * 0.3, ymax * 5.);
            h->GetXaxis()->SetTitleSize(0.05);
            h->GetYaxis()->SetTitleSize(0.05);
            h->GetXaxis()->SetTitleOffset(1.1);
            h->GetYaxis()->SetTitleOffset(1.5);
            if (ikn == 0) h->Draw("E");
            else          h->Draw("E same");
        }

        TLegend* leg1 = new TLegend(0.67, 0.53, 1.08, 0.92);
        leg1->SetBorderSize(0);
        leg1->SetFillStyle(0);
        leg1->SetTextSize(0.034);
        for (int ikn = 0; ikn < nkn; ikn++)
            leg1->AddEntry(hists[ikn], kn_labels[ikn].c_str(), "lep");
        leg1->Draw();

        TLatex lat1;
        lat1.SetNDC();
        lat1.SetTextSize(0.038);
        const std::string label = (ivar == 0) ? "truth p_{T}^{pair}" : "reco p_{T}^{pair}";
        lat1.DrawLatex(0.17, 0.92, ("Pythia fullsim pp24, single-b signal, " + label).c_str());

        // --- Right: stack ---
        c->cd(2);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14);
        gPad->SetLogx();
        gPad->SetLogy();

        THStack* hs = new THStack(("hs_" + outname).c_str(), "");
        for (int ikn = 0; ikn < nkn; ikn++)   // kn0 at bottom
            hs->Add(hists[ikn]);

        hs->Draw("hist");
        hs->GetXaxis()->SetTitle(xtitle.c_str());
        hs->GetYaxis()->SetTitle("d#sigma/dp_{T} [#mub/GeV]");
        hs->GetXaxis()->SetRangeUser(xmin, xmax);
        hs->GetXaxis()->SetTitleSize(0.05);
        hs->GetYaxis()->SetTitleSize(0.05);
        hs->GetXaxis()->SetTitleOffset(1.1);
        hs->GetYaxis()->SetTitleOffset(1.5);

        TLegend* leg2 = new TLegend(0.67, 0.53, 1.08, 0.92);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        leg2->SetTextSize(0.034);
        for (int ikn = 0; ikn < nkn; ikn++)
            leg2->AddEntry(hists[ikn], kn_labels[ikn].c_str(), "f");
        leg2->Draw();

        TLatex lat2;
        lat2.SetNDC();
        lat2.SetTextSize(0.038);
        lat2.DrawLatex(0.17, 0.92, ("Pythia fullsim pp24, single-b signal, " + label).c_str());

        c->SaveAs((output_dir + outname + ".png").c_str());

        for (auto* h : hists) delete h;
        delete hs;
        delete c;
    }
}

// scale_factors[ikn] = N_full / N_test (40k) for each pTHat bin
// err_full = err_test / sqrt(scale_factor)
void plot_stat_error_forecast(int nbins_arg, double xmax_arg, const std::string& suffix,
                               const std::string& subdir = "", double sf_kn5 = 8.) {

    const std::string input_file =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/"
        "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts.root";
    const std::string output_dir =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/plots/";
    gSystem->mkdir((output_dir + subdir).c_str(), true);

    const int nkn = 6;
    const std::array<std::string, nkn> kn_labels = {
        "#hat{p}_{T} 8-14 GeV",
        "#hat{p}_{T} 14-24 GeV",
        "#hat{p}_{T} 24-40 GeV",
        "#hat{p}_{T} 40-70 GeV",
        "#hat{p}_{T} 70-125 GeV",
        "#hat{p}_{T} 125-300 GeV"
    };
    const std::array<int, nkn> colors = {
        kRed+1, kOrange+1, kGreen+2, kCyan+2, kBlue+1, kViolet+1
    };
    const std::array<double, nkn> scale_factors = {51., 90., 60., 30., 8., sf_kn5};

    const int    nbins = nbins_arg;
    const double xmin  = 8., xmax = xmax_arg;
    std::vector<double> edges(nbins + 1);
    const double lmin = std::log(xmin), lmax = std::log(xmax);
    for (int i = 0; i <= nbins; i++)
        edges[i] = std::exp(lmin + i * (lmax - lmin) / nbins);

    const std::array<std::string, 2> var_names   = { "truth_pair_pt", "pair_pt" };
    const std::array<std::string, 2> var_titles  = {
        "truth p_{T}^{pair} [GeV]", "reco p_{T}^{pair} [GeV]"
    };
    const std::array<std::string, 2> filter_strs = {
        "from_same_b",
        "from_same_b && pair_pass_medium"
    };
    const std::array<std::string, 2> label_strs  = {
        "truth p_{T}^{pair}", "reco p_{T}^{pair}"
    };
    const std::array<std::string, 2> out_names   = {
        "truth_pair_pt_kn_stat_error" + suffix,
        "reco_pair_pt_kn_stat_error"  + suffix
    };

    for (int ivar = 0; ivar < 2; ivar++) {
        const auto& var     = var_names[ivar];
        const auto& xtitle  = var_titles[ivar];
        const auto& filter  = filter_strs[ivar];
        const auto& label   = label_strs[ivar];
        const auto& outname = out_names[ivar];

        std::vector<TH1D*> hists(nkn);
        for (int ikn = 0; ikn < nkn; ikn++) {
            const std::string tree = "muon_pair_tree_kin" + std::to_string(ikn) + "_sign2";
            ROOT::RDataFrame df(tree, input_file);
            auto hptr = df.Filter(filter)
                          .Histo1D(ROOT::RDF::TH1DModel{
                              ("hse_" + var + "_kn" + std::to_string(ikn)).c_str(), "",
                              nbins, edges.data()
                          }, var, "weight");
            hists[ikn] = (TH1D*)hptr->Clone();
            hists[ikn]->SetDirectory(nullptr);
            hists[ikn]->Scale(1., "width");
            for (int ib = 1; ib <= hists[ikn]->GetNbinsX(); ib++)
                hists[ikn]->SetBinError(ib, hists[ikn]->GetBinError(ib) / std::sqrt(scale_factors[ikn]));

            hists[ikn]->SetLineColor(colors[ikn]);
            hists[ikn]->SetMarkerColor(colors[ikn]);
            hists[ikn]->SetMarkerStyle(20);
            hists[ikn]->SetMarkerSize(0.7);
            hists[ikn]->SetLineWidth(1);
        }

        // Build relative-error histograms: bin content = err/content = 1/sqrt(N_full_bin)
        std::vector<TH1D*> herr(nkn);
        for (int ikn = 0; ikn < nkn; ikn++) {
            herr[ikn] = (TH1D*)hists[ikn]->Clone(
                ("herr_" + var + "_kn" + std::to_string(ikn) + suffix).c_str());
            herr[ikn]->SetDirectory(nullptr);
            herr[ikn]->Reset();
            for (int ib = 1; ib <= hists[ikn]->GetNbinsX(); ib++) {
                const double c = hists[ikn]->GetBinContent(ib);
                const double e = hists[ikn]->GetBinError(ib);
                herr[ikn]->SetBinContent(ib, c > 0 ? e / c : 0.);
                herr[ikn]->SetBinError(ib, 0.);
            }
            herr[ikn]->SetLineColor(colors[ikn]);
            herr[ikn]->SetLineWidth(2);
        }

        // global y ranges
        double ymax_dist = 0., ymin_dist = 1e30;
        for (auto* h : hists) {
            ymax_dist = std::max(ymax_dist, h->GetMaximum());
            for (int ib = 1; ib <= h->GetNbinsX(); ib++)
                if (h->GetBinContent(ib) > 0)
                    ymin_dist = std::min(ymin_dist, h->GetBinContent(ib));
        }
        if (ymin_dist > 1e29) ymin_dist = 1e-12;

        double ymax_err = 0., ymin_err = 1e30;
        for (auto* h : herr) {
            ymax_err = std::max(ymax_err, h->GetMaximum());
            for (int ib = 1; ib <= h->GetNbinsX(); ib++)
                if (h->GetBinContent(ib) > 0)
                    ymin_err = std::min(ymin_err, h->GetBinContent(ib));
        }
        if (ymin_err > 1e29) ymin_err = 1e-12;

        TCanvas* c = new TCanvas(outname.c_str(), "", 1400, 600);
        c->Divide(2, 1);

        // --- Left: pair pT markers ---
        c->cd(1);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14);
        gPad->SetLogx();
        gPad->SetLogy();
        for (int ikn = 0; ikn < nkn; ikn++) {
            auto* h = hists[ikn];
            h->GetXaxis()->SetTitle(xtitle.c_str());
            h->GetYaxis()->SetTitle("d#sigma/dp_{T} [#mub/GeV]");
            h->GetXaxis()->SetRangeUser(xmin, xmax);
            h->GetYaxis()->SetRangeUser(ymin_dist * 0.3, ymax_dist * 5.);
            h->GetXaxis()->SetTitleSize(0.05);
            h->GetYaxis()->SetTitleSize(0.05);
            h->GetXaxis()->SetTitleOffset(1.1);
            h->GetYaxis()->SetTitleOffset(1.5);
            if (ikn == 0) h->Draw("E");
            else          h->Draw("E same");
        }
        TLegend* leg1 = new TLegend(0.67, 0.53, 1.08, 0.92);
        leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextSize(0.034);
        for (int ikn = 0; ikn < nkn; ikn++)
            leg1->AddEntry(hists[ikn], kn_labels[ikn].c_str(), "lep");
        leg1->Draw();
        TLatex lat1; lat1.SetNDC(); lat1.SetTextSize(0.038);
        lat1.DrawLatex(0.17, 0.92, ("Pythia fullsim pp24, single-b signal, " + label).c_str());

        // --- Right: full-sample stat uncertainty ---
        c->cd(2);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14);
        gPad->SetLogx();
        gPad->SetLogy();
        for (int ikn = 0; ikn < nkn; ikn++) {
            auto* h = herr[ikn];
            h->GetXaxis()->SetTitle(xtitle.c_str());
            h->GetYaxis()->SetTitle("Full-sample rel. stat. error 1/#sqrt{N_{full}}");
            h->GetXaxis()->SetRangeUser(xmin, xmax);
            h->GetYaxis()->SetRangeUser(ymin_err * 0.3, ymax_err * 5.);
            h->GetXaxis()->SetTitleSize(0.05);
            h->GetYaxis()->SetTitleSize(0.05);
            h->GetXaxis()->SetTitleOffset(1.1);
            h->GetYaxis()->SetTitleOffset(1.5);
            if (ikn == 0) h->Draw("hist");
            else          h->Draw("hist same");
        }
        TLegend* leg2 = new TLegend(0.67, 0.53, 1.08, 0.92);
        leg2->SetBorderSize(0); leg2->SetFillStyle(0); leg2->SetTextSize(0.034);
        for (int ikn = 0; ikn < nkn; ikn++)
            leg2->AddEntry(herr[ikn], kn_labels[ikn].c_str(), "l");
        leg2->Draw();
        TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(0.038);
        lat2.DrawLatex(0.17, 0.92, ("Rel. stat. error (full sample), " + label).c_str());

        c->SaveAs((output_dir + subdir + outname + ".png").c_str());

        for (auto* h : hists) delete h;
        for (auto* h : herr)  delete h;
        delete c;
    }
}

// 2D color map: x = pair pT, y = kn region, color = fraction_kn(i) of total quadrature error.
// fraction_kn(i) = σ_err_kn(i) / sqrt(Σ_kn σ_err_kn(i)²),  where σ_err_kn = full-sample abs. error.
// Left: distribution markers.  Right: TH2 color map.
void plot_err_fraction_map(int nbins_arg, double xmax_arg, const std::string& suffix,
                            const std::string& subdir = "", double sf_kn5 = 8.) {

    const std::string input_file =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/"
        "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts.root";
    const std::string output_dir =
        "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/plots/";
    gSystem->mkdir((output_dir + subdir).c_str(), true);

    const int nkn = 6;
    const std::array<std::string, nkn> kn_labels = {
        "#hat{p}_{T} 8-14 GeV",
        "#hat{p}_{T} 14-24 GeV",
        "#hat{p}_{T} 24-40 GeV",
        "#hat{p}_{T} 40-70 GeV",
        "#hat{p}_{T} 70-125 GeV",
        "#hat{p}_{T} 125-300 GeV"
    };
    const std::array<std::string, nkn> kn_short = {
        "8-14", "14-24", "24-40", "40-70", "70-125", "125-300"
    };
    const std::array<int, nkn> colors = {
        kRed+1, kOrange+1, kGreen+2, kCyan+2, kBlue+1, kViolet+1
    };
    const std::array<double, nkn> scale_factors = {51., 90., 60., 30., 8., sf_kn5};

    const int    nbins = nbins_arg;
    const double xmin  = 8., xmax = xmax_arg;
    std::vector<double> edges(nbins + 1);
    const double lmin = std::log(xmin), lmax = std::log(xmax);
    for (int i = 0; i <= nbins; i++)
        edges[i] = std::exp(lmin + i * (lmax - lmin) / nbins);

    const std::array<std::string, 2> var_names   = { "truth_pair_pt", "pair_pt" };
    const std::array<std::string, 2> var_titles  = {
        "truth p_{T}^{pair} [GeV]", "reco p_{T}^{pair} [GeV]"
    };
    const std::array<std::string, 2> filter_strs = {
        "from_same_b",
        "from_same_b && pair_pass_medium"
    };
    const std::array<std::string, 2> label_strs  = {
        "truth p_{T}^{pair}", "reco p_{T}^{pair}"
    };
    const std::array<std::string, 2> out_names   = {
        "truth_pair_pt_kn_err_frac" + suffix,
        "reco_pair_pt_kn_err_frac"  + suffix
    };

    for (int ivar = 0; ivar < 2; ivar++) {
        const auto& var     = var_names[ivar];
        const auto& xtitle  = var_titles[ivar];
        const auto& filter  = filter_strs[ivar];
        const auto& label   = label_strs[ivar];
        const auto& outname = out_names[ivar];

        // Fill differential histograms
        std::vector<TH1D*> hists(nkn);
        for (int ikn = 0; ikn < nkn; ikn++) {
            const std::string tree = "muon_pair_tree_kin" + std::to_string(ikn) + "_sign2";
            ROOT::RDataFrame df(tree, input_file);
            auto hptr = df.Filter(filter)
                          .Histo1D(ROOT::RDF::TH1DModel{
                              ("hef_" + var + "_kn" + std::to_string(ikn)).c_str(), "",
                              nbins, edges.data()
                          }, var, "weight");
            hists[ikn] = (TH1D*)hptr->Clone();
            hists[ikn]->SetDirectory(nullptr);
            hists[ikn]->Scale(1., "width");
            for (int ib = 1; ib <= hists[ikn]->GetNbinsX(); ib++)
                hists[ikn]->SetBinError(ib, hists[ikn]->GetBinError(ib) / std::sqrt(scale_factors[ikn]));

            hists[ikn]->SetLineColor(colors[ikn]);
            hists[ikn]->SetMarkerColor(colors[ikn]);
            hists[ikn]->SetMarkerStyle(20);
            hists[ikn]->SetMarkerSize(0.7);
            hists[ikn]->SetLineWidth(1);
        }

        // Build 2D fraction map: x=pair pT bin, y=kn
        TH2D* h2 = new TH2D(("h2ef_" + outname).c_str(), "",
                             nbins, edges.data(), nkn, -0.5, nkn - 0.5);
        h2->SetDirectory(nullptr);
        for (int ib = 1; ib <= nbins; ib++) {
            double tot2 = 0.;
            for (int ikn = 0; ikn < nkn; ikn++) {
                double e = hists[ikn]->GetBinError(ib);
                tot2 += e * e;
            }
            const double tot = std::sqrt(tot2);
            for (int ikn = 0; ikn < nkn; ikn++) {
                double e = hists[ikn]->GetBinError(ib);
                h2->SetBinContent(ib, ikn + 1, tot > 0 ? e / tot : 0.);
            }
        }
        for (int ikn = 0; ikn < nkn; ikn++)
            h2->GetYaxis()->SetBinLabel(ikn + 1, kn_short[ikn].c_str());

        // y range for left panel
        double ymax_dist = 0., ymin_dist = 1e30;
        for (auto* h : hists) {
            ymax_dist = std::max(ymax_dist, h->GetMaximum());
            for (int ib = 1; ib <= h->GetNbinsX(); ib++)
                if (h->GetBinContent(ib) > 0)
                    ymin_dist = std::min(ymin_dist, h->GetBinContent(ib));
        }
        if (ymin_dist > 1e29) ymin_dist = 1e-12;

        TCanvas* c = new TCanvas(outname.c_str(), "", 1400, 600);
        c->Divide(2, 1);

        // --- Left: distribution markers ---
        c->cd(1);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14);
        gPad->SetLogx();
        gPad->SetLogy();
        for (int ikn = 0; ikn < nkn; ikn++) {
            auto* h = hists[ikn];
            h->GetXaxis()->SetTitle(xtitle.c_str());
            h->GetYaxis()->SetTitle("d#sigma/dp_{T} [#mub/GeV]");
            h->GetXaxis()->SetRangeUser(xmin, xmax);
            h->GetYaxis()->SetRangeUser(ymin_dist * 0.3, ymax_dist * 5.);
            h->GetXaxis()->SetTitleSize(0.05);
            h->GetYaxis()->SetTitleSize(0.05);
            h->GetXaxis()->SetTitleOffset(1.1);
            h->GetYaxis()->SetTitleOffset(1.5);
            if (ikn == 0) h->Draw("E");
            else          h->Draw("E same");
        }
        TLegend* leg1 = new TLegend(0.67, 0.53, 1.08, 0.92);
        leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextSize(0.034);
        for (int ikn = 0; ikn < nkn; ikn++)
            leg1->AddEntry(hists[ikn], kn_labels[ikn].c_str(), "lep");
        leg1->Draw();
        TLatex lat1; lat1.SetNDC(); lat1.SetTextSize(0.038);
        lat1.DrawLatex(0.17, 0.92, ("Pythia fullsim pp24, single-b signal, " + label).c_str());

        // --- Right: fraction map ---
        c->cd(2);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.16);   // room for z-axis color bar
        gPad->SetBottomMargin(0.14);
        gPad->SetLogx();
        gStyle->SetPalette(kViridis);
        h2->GetXaxis()->SetTitle(xtitle.c_str());
        h2->GetYaxis()->SetTitle("#hat{p}_{T} range [GeV]");
        h2->GetZaxis()->SetTitle("#sigma_{err,kn} / #sigma_{err,total}");
        h2->GetXaxis()->SetRangeUser(xmin, xmax);
        h2->SetMinimum(0.);
        h2->SetMaximum(1.);
        h2->GetXaxis()->SetTitleSize(0.05);
        h2->GetYaxis()->SetTitleSize(0.05);
        h2->GetZaxis()->SetTitleSize(0.045);
        h2->GetXaxis()->SetTitleOffset(1.1);
        h2->GetYaxis()->SetTitleOffset(1.2);
        h2->GetZaxis()->SetTitleOffset(1.3);
        h2->GetYaxis()->SetLabelSize(0.05);
        h2->Draw("colz");
        TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(0.038);
        lat2.DrawLatex(0.15, 0.92, ("Full-sample error fraction, " + label).c_str());

        c->SaveAs((output_dir + subdir + outname + ".png").c_str());

        for (auto* h : hists) delete h;
        delete h2;
        delete c;
    }
}

void replot_scale_forecast() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // new scales (kn5 sf=8)
    plot_stat_error_forecast(20, 120., "");
    plot_stat_error_forecast(25, 150., "_150GeV");
    plot_err_fraction_map(20, 120., "");
    plot_err_fraction_map(25, 150., "_150GeV");
    // old scales (kn5 sf=1) → old_scales/
    plot_stat_error_forecast(20, 120., "", "old_scales/", 1.);
    plot_stat_error_forecast(25, 150., "_150GeV", "old_scales/", 1.);
    plot_err_fraction_map(20, 120., "", "old_scales/", 1.);
    plot_err_fraction_map(25, 150., "_150GeV", "old_scales/", 1.);
}

void plot_pythia_fullsim_kn_pt_crossx() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    plot_impl(20, 120., "");
    plot_impl(25, 150., "_150GeV");
    plot_stat_error_forecast(20, 120., "");
    plot_stat_error_forecast(25, 150., "_150GeV");
    plot_err_fraction_map(20, 120., "");
    plot_err_fraction_map(25, 150., "_150GeV");
}
