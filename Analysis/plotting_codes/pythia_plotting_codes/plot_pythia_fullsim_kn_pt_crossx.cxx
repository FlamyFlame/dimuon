// Differential cross-section vs pair pT per pTHat (kn) range.
// Selection: muon_pair_tree_kin*_sign2 with from_same_b (nominal), or -- ss_loose_mass=true --
//            muon_pair_tree_kin*_sign1 with truth/reco m_mumu < MCRequest::kLooseMassMax.
// Two plots: truth_pair_pt and reco pair_pt (pair_pass_medium additionally required).
// Left subplot: markers+errorbars per kn.  Right subplot: stack.
// Entry point calls two binning versions: 20 bins 9-120 GeV and 25 bins 9-150 GeV.
// (The low edge follows ParamsSet::signal_pair_pt_min, 8 -> 9 on 2026-09-08.)

// UNITS (do not "fix" this back): the per-pair `weight` is sigma*genFiltEff*r_isospin/N_beam
// with the AMI cross-section in **nb** (Analysis/docs/ami_weights.md). No unit conversion is
// applied anywhere below, so dsigma/dpT is in **nb/GeV**. The axis was previously labelled
// [ub/GeV] -- a 1000x mislabel. Comparing to pp DATA (dsigma = N/L, L in pb^-1) needs nb->pb, x1000.

#include "../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../Utilities/PtHatKn45ProjectedStats.h"
#include "../../Utilities/MCRequestLooseMassFilter.h"
#include <TString.h>   // Form(), used by the population helpers below
#include <string>

// ============================ CONFIG ============================
// g_is_test_sample : which pp24-fullsim production to read.
//   true  = TEST sample -- 4 isospin beams (produced that way BY MISTAKE), so its absolute sigma
//           carries the Pb 4:6:6:9 isospin AVERAGE and is NOT a physical pp cross-section.
//   false = FULL sample -- pp beam only, isospin weight 1 => an HONEST pp cross-section.
//           Files end in "_full"; plots go to that sample's own plots/ dir (no clobber).
//   One switch -- it also fixes the isospin treatment upstream (FullSimSampleType.h).
// g_use_tight_wp : NOMINAL muon WP is TIGHT (feedback_plots_wp_config_var, muon_wp_registry.md).
//   This macro previously HARD-CODED Medium while every other pp-fullsim stage used Tight.
// NOT `static`: ACLiC-compiled internal-linkage globals are invisible to the ROOT interpreter,
// so a caller could not override them and every run silently used the defaults. The entry point
// below takes them as arguments and sets these.
// g_ss_loose_mass : which PAIR POPULATION the per-slice statistics plots describe.
//   false = NOMINAL: single-b signal pairs -- opposite-sign tree, `from_same_b` (truth-level, so
//           the same filter serves the truth and the reco variable).
//   true  = SAME-SIGN pairs under the LOOSE MASS cut of the additional-statistics MC request
//           (docs/tracking/pthat_slice_mass_stats_sample_request.md D9, 2026-09-17): same-sign
//           tree (`_sign1`, split by truth_same_sign upstream), truth pair pT under
//           truth_minv < kLooseMassMax and reco pair pT under minv < kLooseMassMax -- each variable
//           cut at ITS OWN level, because the request's generator filter is a truth-mass cut while
//           the analysis consumes reco mass. Outputs go to their own subdirectory (ss_mass_lt10/)
//           with the nominal file names; the nominal PNGs are never touched. The projected
//           function is nominal-only and throws on this flag.
bool g_is_test_sample = true;
bool g_use_tight_wp   = true;
bool g_ss_loose_mass  = false;

static std::string SampleSuffix() { return g_is_test_sample ? "" : "_full"; }
static std::string SampleDir()    { return FullSimSampleInputDir(FullSimSampleType::pp, g_is_test_sample); }
static std::string InputFile()    { return SampleDir() + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts"
                                           + SampleSuffix() + ".root"; }
static std::string OutputDir()    { return SampleDir() + "plots/" + (g_ss_loose_mass ? "ss_mass_lt10/" : ""); }
static std::string PairWPFilter() { return g_use_tight_wp ? "pair_pass_tight" : "pair_pass_medium"; }
// The pair population (see g_ss_loose_mass): tree, per-level filters, and the label naming it.
static std::string PairTree(int ikn) { return "muon_pair_tree_kin" + std::to_string(ikn)
                                              + (g_ss_loose_mass ? "_sign1" : "_sign2"); }
static std::string TruthFilter()     { return g_ss_loose_mass
                                              ? std::string(Form("truth_minv < %g", MCRequest::kLooseMassMax))
                                              : std::string("from_same_b"); }
static std::string RecoFilter()      { return (g_ss_loose_mass
                                               ? std::string(Form("minv < %g", MCRequest::kLooseMassMax))
                                               : std::string("from_same_b")) + " && " + PairWPFilter(); }
// The same-sign header is ~14 characters longer than the nominal one, which fills the panel
// width exactly at 0.038; it is shrunk so the header cannot run under the legend.
static double      HeaderTextSize()  { return g_ss_loose_mass ? 0.032 : 0.038; }
static std::string SelLabel()        { return g_ss_loose_mass
                                              ? std::string(Form("same sign, m_{#mu#mu} < %g GeV", MCRequest::kLooseMassMax))
                                              : std::string("single-b signal"); }

// STATISTICAL-FORECAST factors N_full/N_test. Meaningful ONLY on the TEST sample, where they answer
// "what will the error bar be once the full sample exists" (they are applied to SetBinError ONLY;
// central values are untouched). On the FULL sample the histogram errors ALREADY reflect the real
// statistics, so dividing by sqrt(sf) again would understate them a SECOND time -- by a further
// 7.1x / 9.5x / 7.7x / 5.5x / 2.8x, SLICE-DEPENDENTLY (cancels nowhere) -- while the cross-section
// curve still looked perfect. A SILENT error. Hence the forecast is a NO-OP on the full sample.
static double ForecastScale(double sf_test) { return g_is_test_sample ? sf_test : 1.0; }
// ================================================================

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
#include <stdexcept>
#include <string>
#include <vector>
#include <array>
#include "../../MuonObjectsParamsAndHelpers/ParamsSet.h"   // signal_pair_pt_min / pTbins -- axes read, never retyped

void plot_impl(int nbins_arg, double xmax_arg, const std::string& suffix) {

    const std::string input_file = InputFile();
    const std::string output_dir = OutputDir();
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
    const double xmin  = ParamsSet::signal_pair_pt_min, xmax = xmax_arg;  // 8 -> 9 (2026-09-08)
    std::vector<double> edges(nbins + 1);
    const double lmin = std::log(xmin), lmax = std::log(xmax);
    for (int i = 0; i <= nbins; i++)
        edges[i] = std::exp(lmin + i * (lmax - lmin) / nbins);

    // Two variables: [0] truth_pair_pt, [1] pair_pt (reco)
    const std::array<std::string, 2> var_names   = { "truth_pair_pt", "pair_pt" };
    const std::array<std::string, 2> var_titles  = {
        "truth p_{T}^{pair} [GeV]", "reco p_{T}^{pair} [GeV]"
    };
    const std::array<std::string, 2> filter_strs = { TruthFilter(), RecoFilter() };
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
            ROOT::RDataFrame df(PairTree(ikn), input_file);
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
            h->GetYaxis()->SetTitle("d#sigma/dp_{T} [nb/GeV]");
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
        lat1.SetTextSize(HeaderTextSize());
        const std::string label = (ivar == 0) ? "truth p_{T}^{pair}" : "reco p_{T}^{pair}";
        // HONESTY (key physics observable). The TEST sample carries the Pb 4:6:6:9 isospin AVERAGE
        // (it was produced with 4 beams by mistake), so its absolute sigma is NOT a physical pp
        // cross-section -- say so. The FULL sample is pp-beam-only with isospin weight 1, so its
        // sigma IS an honest pp cross-section and the warning must NOT be printed (it would be a
        // FALSE warning).
        lat1.DrawLatex(0.17, 0.92, (std::string("Pythia fullsim pp24 ")
            + (g_is_test_sample ? "TEST sample" : "FULL sample")
            + ", " + SelLabel() + ", " + label).c_str());
        if (g_is_test_sample) {
            TLatex lat_warn;
            lat_warn.SetNDC();
            lat_warn.SetTextSize(0.026);
            lat_warn.SetTextColor(kRed + 1);
            // Bottom-left: the spectrum falls away from this corner, so it clears both the curves
            // and the legend (at 0.92 it collided with the legend and was clipped).
            lat_warn.DrawLatex(0.20, 0.235, "Pb isospin avg (4:6:6:9) applied to a pp sample");
            lat_warn.DrawLatex(0.20, 0.200, "#Rightarrow NOT a physical pp #sigma");
        }

        // THStack on LOG y -- the WIDE-DYNAMIC-RANGE EXCEPTION (.claude/conventions/atlas-plotting.md
        // "Stacked histograms"; criterion C7 of physics-results-review.md). A THStack defaults to a
        // LINEAR y-axis, because a band's THICKNESS should be proportional to its contribution.
        // Deliberately overridden here: dsigma/dpT per pT-hat slice spans ~5 DECADES, so on a linear axis
        // every slice but the lowest collapses onto zero and the whole high-pT tail -- the physics of
        // interest -- becomes an invisible flat line. Unreadable is worse than distorted.
        // Consequence accepted: on log y the band thickness is NOT proportional to the contribution, so
        // this stack must NOT be read as "parts of a whole"/fractions -- use the left (marker) panel for
        // per-slice values. The magnitude-ordering rule is axis-independent and still applies.
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
        hs->GetYaxis()->SetTitle("d#sigma/dp_{T} [nb/GeV]");
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
        lat2.SetTextSize(HeaderTextSize());
        lat2.DrawLatex(0.17, 0.92, (std::string("Pythia fullsim pp24 ")
            + (g_is_test_sample ? "TEST sample" : "FULL sample")
            + ", " + SelLabel() + ", " + label).c_str());
        if (g_is_test_sample) {
            TLatex lat2_warn;
            lat2_warn.SetNDC();
            lat2_warn.SetTextSize(0.026);
            lat2_warn.SetTextColor(kRed + 1);
            // On the STACK panel the bands fill the bottom-left, so the warning goes just under the
            // title, left of the legend, where the stack has not yet risen.
            lat2_warn.DrawLatex(0.20, 0.865, "Pb isospin avg (4:6:6:9) applied to a pp sample");
            lat2_warn.DrawLatex(0.20, 0.830, "#Rightarrow NOT a physical pp #sigma");
        }

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

    const std::string input_file = InputFile();
    const std::string output_dir = OutputDir();
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
    const double xmin  = ParamsSet::signal_pair_pt_min, xmax = xmax_arg;  // 8 -> 9 (2026-09-08)
    std::vector<double> edges(nbins + 1);
    const double lmin = std::log(xmin), lmax = std::log(xmax);
    for (int i = 0; i <= nbins; i++)
        edges[i] = std::exp(lmin + i * (lmax - lmin) / nbins);

    const std::array<std::string, 2> var_names   = { "truth_pair_pt", "pair_pt" };
    const std::array<std::string, 2> var_titles  = {
        "truth p_{T}^{pair} [GeV]", "reco p_{T}^{pair} [GeV]"
    };
    const std::array<std::string, 2> filter_strs = { TruthFilter(), RecoFilter() };
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
            ROOT::RDataFrame df(PairTree(ikn), input_file);
            auto hptr = df.Filter(filter)
                          .Histo1D(ROOT::RDF::TH1DModel{
                              ("hse_" + var + "_kn" + std::to_string(ikn)).c_str(), "",
                              nbins, edges.data()
                          }, var, "weight");
            hists[ikn] = (TH1D*)hptr->Clone();
            hists[ikn]->SetDirectory(nullptr);
            hists[ikn]->Scale(1., "width");
            for (int ib = 1; ib <= hists[ikn]->GetNbinsX(); ib++)
                hists[ikn]->SetBinError(ib, hists[ikn]->GetBinError(ib) / std::sqrt(ForecastScale(scale_factors[ikn])));

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
            h->GetYaxis()->SetTitle("d#sigma/dp_{T} [nb/GeV]");
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
        TLatex lat1; lat1.SetNDC(); lat1.SetTextSize(HeaderTextSize());
        // HONESTY (key physics observable). The TEST sample carries the Pb 4:6:6:9 isospin AVERAGE
        // (it was produced with 4 beams by mistake), so its absolute sigma is NOT a physical pp
        // cross-section -- say so. The FULL sample is pp-beam-only with isospin weight 1, so its
        // sigma IS an honest pp cross-section and the warning must NOT be printed (it would be a
        // FALSE warning).
        lat1.DrawLatex(0.17, 0.92, (std::string("Pythia fullsim pp24 ")
            + (g_is_test_sample ? "TEST sample" : "FULL sample")
            + ", " + SelLabel() + ", " + label).c_str());
        if (g_is_test_sample) {
            TLatex lat_warn;
            lat_warn.SetNDC();
            lat_warn.SetTextSize(0.026);
            lat_warn.SetTextColor(kRed + 1);
            // Bottom-left: the spectrum falls away from this corner, so it clears both the curves
            // and the legend (at 0.92 it collided with the legend and was clipped).
            lat_warn.DrawLatex(0.20, 0.235, "Pb isospin avg (4:6:6:9) applied to a pp sample");
            lat_warn.DrawLatex(0.20, 0.200, "#Rightarrow NOT a physical pp #sigma");
        }

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
        TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(HeaderTextSize());
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

    const std::string input_file = InputFile();
    const std::string output_dir = OutputDir();
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
    const double xmin  = ParamsSet::signal_pair_pt_min, xmax = xmax_arg;  // 8 -> 9 (2026-09-08)
    std::vector<double> edges(nbins + 1);
    const double lmin = std::log(xmin), lmax = std::log(xmax);
    for (int i = 0; i <= nbins; i++)
        edges[i] = std::exp(lmin + i * (lmax - lmin) / nbins);

    const std::array<std::string, 2> var_names   = { "truth_pair_pt", "pair_pt" };
    const std::array<std::string, 2> var_titles  = {
        "truth p_{T}^{pair} [GeV]", "reco p_{T}^{pair} [GeV]"
    };
    const std::array<std::string, 2> filter_strs = { TruthFilter(), RecoFilter() };
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
            ROOT::RDataFrame df(PairTree(ikn), input_file);
            auto hptr = df.Filter(filter)
                          .Histo1D(ROOT::RDF::TH1DModel{
                              ("hef_" + var + "_kn" + std::to_string(ikn)).c_str(), "",
                              nbins, edges.data()
                          }, var, "weight");
            hists[ikn] = (TH1D*)hptr->Clone();
            hists[ikn]->SetDirectory(nullptr);
            hists[ikn]->Scale(1., "width");
            for (int ib = 1; ib <= hists[ikn]->GetNbinsX(); ib++)
                hists[ikn]->SetBinError(ib, hists[ikn]->GetBinError(ib) / std::sqrt(ForecastScale(scale_factors[ikn])));

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
            h->GetYaxis()->SetTitle("d#sigma/dp_{T} [nb/GeV]");
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
        TLatex lat1; lat1.SetNDC(); lat1.SetTextSize(HeaderTextSize());
        // HONESTY (key physics observable). The TEST sample carries the Pb 4:6:6:9 isospin AVERAGE
        // (it was produced with 4 beams by mistake), so its absolute sigma is NOT a physical pp
        // cross-section -- say so. The FULL sample is pp-beam-only with isospin weight 1, so its
        // sigma IS an honest pp cross-section and the warning must NOT be printed (it would be a
        // FALSE warning).
        lat1.DrawLatex(0.17, 0.92, (std::string("Pythia fullsim pp24 ")
            + (g_is_test_sample ? "TEST sample" : "FULL sample")
            + ", " + SelLabel() + ", " + label).c_str());
        if (g_is_test_sample) {
            TLatex lat_warn;
            lat_warn.SetNDC();
            lat_warn.SetTextSize(0.026);
            lat_warn.SetTextColor(kRed + 1);
            // Bottom-left: the spectrum falls away from this corner, so it clears both the curves
            // and the legend (at 0.92 it collided with the legend and was clipped).
            lat_warn.DrawLatex(0.20, 0.235, "Pb isospin avg (4:6:6:9) applied to a pp sample");
            lat_warn.DrawLatex(0.20, 0.200, "#Rightarrow NOT a physical pp #sigma");
        }

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
        TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(HeaderTextSize());
        lat2.DrawLatex(0.15, 0.92, ("Full-sample error fraction, " + label).c_str());

        c->SaveAs((output_dir + subdir + outname + ".png").c_str());

        for (auto* h : hists) delete h;
        delete h2;
        delete c;
    }
}

// Double-ratio color map: (err fraction) / (xsec fraction) per (pair pT bin, kn).
// Value > 1: kn contributes more to error than to statistics → statistical bottleneck.
// Value < 1: kn is over-sampled relative to its cross-section share.
void plot_err_ratio_map(int nbins_arg, double xmax_arg, const std::string& suffix,
                         const std::string& subdir = "", double sf_kn5 = 8.) {

    const std::string input_file = InputFile();
    const std::string output_dir = OutputDir();
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
    const double xmin  = ParamsSet::signal_pair_pt_min, xmax = xmax_arg;  // 8 -> 9 (2026-09-08)
    std::vector<double> edges(nbins + 1);
    const double lmin = std::log(xmin), lmax = std::log(xmax);
    for (int i = 0; i <= nbins; i++)
        edges[i] = std::exp(lmin + i * (lmax - lmin) / nbins);

    const std::array<std::string, 2> var_names   = { "truth_pair_pt", "pair_pt" };
    const std::array<std::string, 2> var_titles  = {
        "truth p_{T}^{pair} [GeV]", "reco p_{T}^{pair} [GeV]"
    };
    const std::array<std::string, 2> filter_strs = { TruthFilter(), RecoFilter() };
    const std::array<std::string, 2> label_strs  = {
        "truth p_{T}^{pair}", "reco p_{T}^{pair}"
    };
    const std::array<std::string, 2> out_names   = {
        "truth_pair_pt_kn_err_ratio" + suffix,
        "reco_pair_pt_kn_err_ratio"  + suffix
    };

    for (int ivar = 0; ivar < 2; ivar++) {
        const auto& var     = var_names[ivar];
        const auto& xtitle  = var_titles[ivar];
        const auto& filter  = filter_strs[ivar];
        const auto& label   = label_strs[ivar];
        const auto& outname = out_names[ivar];

        // Fill differential histograms with full-sample error bars
        std::vector<TH1D*> hists(nkn);
        for (int ikn = 0; ikn < nkn; ikn++) {
            ROOT::RDataFrame df(PairTree(ikn), input_file);
            auto hptr = df.Filter(filter)
                          .Histo1D(ROOT::RDF::TH1DModel{
                              ("her_" + var + "_kn" + std::to_string(ikn)).c_str(), "",
                              nbins, edges.data()
                          }, var, "weight");
            hists[ikn] = (TH1D*)hptr->Clone();
            hists[ikn]->SetDirectory(nullptr);
            hists[ikn]->Scale(1., "width");
            for (int ib = 1; ib <= hists[ikn]->GetNbinsX(); ib++)
                hists[ikn]->SetBinError(ib, hists[ikn]->GetBinError(ib) / std::sqrt(ForecastScale(scale_factors[ikn])));

            hists[ikn]->SetLineColor(colors[ikn]);
            hists[ikn]->SetMarkerColor(colors[ikn]);
            hists[ikn]->SetMarkerStyle(20);
            hists[ikn]->SetMarkerSize(0.7);
            hists[ikn]->SetLineWidth(1);
        }

        // Build 2D double-ratio map: (err_frac_kn) / (xsec_frac_kn)
        TH2D* h2 = new TH2D(("h2er_" + outname).c_str(), "",
                             nbins, edges.data(), nkn, -0.5, nkn - 0.5);
        h2->SetDirectory(nullptr);
        for (int ib = 1; ib <= nbins; ib++) {
            double tot_err2 = 0., tot_xsec = 0.;
            for (int ikn = 0; ikn < nkn; ikn++) {
                double e = hists[ikn]->GetBinError(ib);
                tot_err2 += e * e;
                tot_xsec += hists[ikn]->GetBinContent(ib);
            }
            const double tot_err = std::sqrt(tot_err2);
            for (int ikn = 0; ikn < nkn; ikn++) {
                double e        = hists[ikn]->GetBinError(ib);
                double c        = hists[ikn]->GetBinContent(ib);
                double ferr     = (tot_err  > 0) ? e / tot_err  : 0.;
                double fxsec    = (tot_xsec > 0) ? c / tot_xsec : 0.;
                double ratio    = (fxsec    > 0) ? ferr / fxsec  : 0.;
                h2->SetBinContent(ib, ikn + 1, ratio);
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

        // --- Left: distribution markers with full-sample error bars ---
        c->cd(1);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14);
        gPad->SetLogx();
        gPad->SetLogy();
        for (int ikn = 0; ikn < nkn; ikn++) {
            auto* h = hists[ikn];
            h->GetXaxis()->SetTitle(xtitle.c_str());
            h->GetYaxis()->SetTitle("d#sigma/dp_{T} [nb/GeV]");
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
        TLatex lat1; lat1.SetNDC(); lat1.SetTextSize(HeaderTextSize());
        // HONESTY (key physics observable). The TEST sample carries the Pb 4:6:6:9 isospin AVERAGE
        // (it was produced with 4 beams by mistake), so its absolute sigma is NOT a physical pp
        // cross-section -- say so. The FULL sample is pp-beam-only with isospin weight 1, so its
        // sigma IS an honest pp cross-section and the warning must NOT be printed (it would be a
        // FALSE warning).
        lat1.DrawLatex(0.17, 0.92, (std::string("Pythia fullsim pp24 ")
            + (g_is_test_sample ? "TEST sample" : "FULL sample")
            + ", " + SelLabel() + ", " + label).c_str());
        if (g_is_test_sample) {
            TLatex lat_warn;
            lat_warn.SetNDC();
            lat_warn.SetTextSize(0.026);
            lat_warn.SetTextColor(kRed + 1);
            // Bottom-left: the spectrum falls away from this corner, so it clears both the curves
            // and the legend (at 0.92 it collided with the legend and was clipped).
            lat_warn.DrawLatex(0.20, 0.235, "Pb isospin avg (4:6:6:9) applied to a pp sample");
            lat_warn.DrawLatex(0.20, 0.200, "#Rightarrow NOT a physical pp #sigma");
        }

        // --- Right: double-ratio color map ---
        c->cd(2);
        gPad->SetLeftMargin(0.14);
        gPad->SetRightMargin(0.16);
        gPad->SetBottomMargin(0.14);
        gPad->SetLogx();
        gStyle->SetPalette(kViridis);
        h2->GetXaxis()->SetTitle(xtitle.c_str());
        h2->GetYaxis()->SetTitle("#hat{p}_{T} range [GeV]");
        h2->GetZaxis()->SetTitle("(#sigma_{err,kn}/#sigma_{err,tot}) / (#sigma_{kn}/#sigma_{tot})");
        h2->GetXaxis()->SetRangeUser(xmin, xmax);
        h2->SetMinimum(0.);
        h2->GetXaxis()->SetTitleSize(0.05);
        h2->GetYaxis()->SetTitleSize(0.05);
        h2->GetZaxis()->SetTitleSize(0.038);
        h2->GetXaxis()->SetTitleOffset(1.1);
        h2->GetYaxis()->SetTitleOffset(1.2);
        h2->GetZaxis()->SetTitleOffset(1.5);
        h2->GetYaxis()->SetLabelSize(0.05);
        h2->Draw("colz");
        TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(HeaderTextSize());
        lat2.DrawLatex(0.15, 0.92, ("Full-sample err/xsec ratio, " + label).c_str());

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
    plot_err_ratio_map(20, 120., "");
    plot_err_ratio_map(25, 150., "_150GeV");
    // old scales (kn5 sf=1) → old_scales/
    plot_stat_error_forecast(20, 120., "", "old_scales/", 1.);
    plot_stat_error_forecast(25, 150., "_150GeV", "old_scales/", 1.);
    plot_err_fraction_map(20, 120., "", "old_scales/", 1.);
    plot_err_fraction_map(25, 150., "_150GeV", "old_scales/", 1.);
    plot_err_ratio_map(20, 120., "", "old_scales/", 1.);
    plot_err_ratio_map(25, 150., "_150GeV", "old_scales/", 1.);
}

// Args (NOT interpreter-global assignment — see the g_* note): is_test_sample selects TEST vs
// FULL production (input dir, "_full" suffix, honesty caption, forecast-factor no-op);
// use_tight_wp selects the pair WP (nominal TIGHT).
void plot_pythia_fullsim_kn_pt_crossx(bool is_test_sample = true, bool use_tight_wp = true,
                                      bool ss_loose_mass = false) {
    g_is_test_sample = is_test_sample;
    g_use_tight_wp   = use_tight_wp;
    g_ss_loose_mass  = ss_loose_mass;
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    plot_impl(20, 120., "");
    plot_impl(25, 150., "_150GeV");
    plot_stat_error_forecast(20, 120., "");
    plot_stat_error_forecast(25, 150., "_150GeV");
    plot_err_fraction_map(20, 120., "");
    plot_err_fraction_map(25, 150., "_150GeV");
}

// ================================================================================================
// PROJECTED-STATISTICS companion plot (2026-09-03, user request; see
// docs/tracking/pythia_pp24_pthat_stats_projection.md). Answers "what would this spectrum look
// like if the two highest pT-hat slices (kn4 = 70-125 GeV, kn5 = 125-300 GeV) had N_TARGET events
// instead of the N_CURRENT the FULL sample was actually produced with" -- a decision aid for an
// MC production request, NOT a new measurement.
//
// PHYSICS (doc Physics Procedure §1): the per-pair weight w = sigma_slice*genFiltEff/N_slice
// scales as 1/N while the pair COUNT scales as N, so sum(w) -- the bin content plotted here -- is
// an UNBIASED estimator of the slice's cross section independent of N. Only its STATISTICAL
// UNCERTAINTY shrinks, as 1/sqrt(N). So: central values (bin content) of kn4/kn5 are left
// UNTOUCHED, and their bin ERROR is divided by sqrt(sf), sf = N_target/N_current. kn0-3 (not part
// of the request) get sf = 1, i.e. are byte-identical to the nominal `plot_impl` output. This is
// the SAME mechanism as `ForecastScale`/`plot_stat_error_forecast` above (TEST->FULL forecast);
// here it runs FULL-sample -> a LARGER FULL sample, unconditionally (not gated by
// g_is_test_sample), for kn4/kn5 only.
//
// N_current/N_target (measured, both slices alike) live in Utilities/PtHatKn45ProjectedStats.h --
// the SINGLE SOURCE OF TRUTH shared with the projected SS/OS statistics CSVs, so the plot and the
// tables can never disagree about the scale factor.
//
// The ORIGINAL plot_impl PNGs are never touched: this writes to its OWN subdirectory with its own
// file names (user instruction: "keep that original plot, but make a new plot").
void plot_impl_projected(int nbins_arg, double xmax_arg, const std::string& suffix) {

    if (g_is_test_sample)
        throw std::runtime_error("plot_impl_projected: projecting additional statistics is only "
            "meaningful starting from the CURRENT FULL sample (g_is_test_sample must be false)");
    if (g_ss_loose_mass)
        throw std::runtime_error("plot_impl_projected: the projection is defined for the single-b "
            "signal population only (it hard-codes the opposite-sign tree and from_same_b)");

    const double kSf = PtHatKn45Projected::kSf;   // = 3.75001..., both slices alike

    const std::string input_file = InputFile();
    const std::string output_dir = OutputDir() + "projected_stats/";
    gSystem->mkdir(output_dir.c_str(), true);

    const int nkn = 6;
    const std::array<std::string, nkn> kn_labels_nominal = {
        "#hat{p}_{T} 8-14 GeV",
        "#hat{p}_{T} 14-24 GeV",
        "#hat{p}_{T} 24-40 GeV",
        "#hat{p}_{T} 40-70 GeV",
        "#hat{p}_{T} 70-125 GeV",
        "#hat{p}_{T} 125-300 GeV"
    };
    std::array<std::string, nkn> kn_labels = kn_labels_nominal;
    kn_labels[4] += " (proj. 1.2M evt)";
    kn_labels[5] += " (proj. 1.2M evt)";
    const std::array<int, nkn> colors = {
        kRed+1, kOrange+1, kGreen+2, kCyan+2, kBlue+1, kViolet+1
    };

    const int    nbins = nbins_arg;
    const double xmin  = ParamsSet::signal_pair_pt_min, xmax = xmax_arg;  // 8 -> 9 (2026-09-08)
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
        ("from_same_b && " + PairWPFilter())
    };
    const std::array<std::string, 2> label_strs  = {
        "truth p_{T}^{pair}", "reco p_{T}^{pair}"
    };
    const std::array<std::string, 2> out_names   = {
        "truth_pair_pt_kn_projected" + suffix, "reco_pair_pt_kn_projected" + suffix
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
                              ("hproj_" + var + "_kn" + std::to_string(ikn)).c_str(), "",
                              nbins, edges.data()
                          }, var, "weight");
            hists[ikn] = (TH1D*)hptr->Clone();
            hists[ikn]->SetDirectory(nullptr);
            hists[ikn]->Scale(1., "width");   // differential; central values UNCHANGED by the projection

            if (ikn == 4 || ikn == 5)         // kn4/kn5 ONLY: stat. error -> what N_target buys
                for (int ib = 1; ib <= hists[ikn]->GetNbinsX(); ib++)
                    hists[ikn]->SetBinError(ib, hists[ikn]->GetBinError(ib) / std::sqrt(kSf));

            hists[ikn]->SetLineColor(colors[ikn]);
            hists[ikn]->SetMarkerColor(colors[ikn]);
            hists[ikn]->SetMarkerStyle(20);
            hists[ikn]->SetMarkerSize(0.7);
            hists[ikn]->SetLineWidth(1);
        }

        double ymax = 0.;
        for (auto* h : hists) ymax = std::max(ymax, h->GetMaximum());
        double ymin_nonzero = 1e30;
        for (auto* h : hists)
            for (int ib = 1; ib <= h->GetNbinsX(); ib++)
                if (h->GetBinContent(ib) > 0)
                    ymin_nonzero = std::min(ymin_nonzero, h->GetBinContent(ib));
        if (ymin_nonzero > 1e29) ymin_nonzero = 1e-12;

        TCanvas* c = new TCanvas(outname.c_str(), "", 800, 700);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.14);
        gPad->SetLogx();
        gPad->SetLogy();

        for (int ikn = 0; ikn < nkn; ikn++) {
            auto* h = hists[ikn];
            h->GetXaxis()->SetTitle(xtitle.c_str());
            h->GetYaxis()->SetTitle("d#sigma/dp_{T} [nb/GeV]");
            h->GetXaxis()->SetRangeUser(xmin, xmax);
            h->GetYaxis()->SetRangeUser(ymin_nonzero * 0.3, ymax * 5.);
            h->GetXaxis()->SetTitleSize(0.045);
            h->GetYaxis()->SetTitleSize(0.045);
            h->GetXaxis()->SetTitleOffset(1.1);
            h->GetYaxis()->SetTitleOffset(1.6);
            if (ikn == 0) h->Draw("E");
            else          h->Draw("E same");
        }

        TLegend* leg = new TLegend(0.50, 0.60, 0.93, 0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.028);
        for (int ikn = 0; ikn < nkn; ikn++)
            leg->AddEntry(hists[ikn], kn_labels[ikn].c_str(), "lep");
        leg->Draw();

        TLatex lat;
        lat.SetNDC();
        lat.SetTextSize(0.032);
        lat.DrawLatex(0.17, 0.93, (std::string("Pythia fullsim pp24 FULL sample, single-b signal, ")
            + label).c_str());
        TLatex lat_note;
        lat_note.SetNDC();
        lat_note.SetTextSize(0.024);
        lat_note.SetTextColor(kBlue+2);
        lat_note.DrawLatex(0.17, 0.895,
            "Projection: kn4/kn5 error bars scaled to 1.2M evt (from 320k); central values unchanged");

        c->SaveAs((output_dir + outname + ".png").c_str());

        for (auto* h : hists) delete h;
        delete c;
    }
}

// Entry point: nominal-WP FULL-sample projection, both binnings. is_test_sample must already be
// false when this runs (set by a prior plot_pythia_fullsim_kn_pt_crossx(false, ...) call, or here
// directly) -- projecting from the TEST sample would answer a different, uninteresting question.
void plot_pythia_fullsim_kn_pt_crossx_projected(bool use_tight_wp = true) {
    g_is_test_sample = false;
    g_use_tight_wp   = use_tight_wp;
    g_ss_loose_mass  = false;
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    plot_impl_projected(20, 120., "");
    plot_impl_projected(25, 150., "_150GeV");
}
