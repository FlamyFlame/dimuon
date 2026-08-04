// =============================================================================
// FitMCSinglesEffcy.cxx
//
// MC single-muon mu4 turn-on fits (docs/tracking/mc_trigger_efficiency.md §3.1),
// consumed by FillMCTrigEffHists Step 3 (§3.3 inverse weights).
//
// Mirrors the DATA fitter Analysis/SingleMuEffcyPtTurnOnFitter.cxx EXACTLY in the
// fit ingredients (same functional forms, parameter init/limits, "QR" fit option,
// [4,60] range, TFormula-string TF1 so the analytic formula persists on Write):
//   pp      -> erf_plus_log  (data pp nominal;   fitter lines 183-219, entry 451-457)
//   overlay -> fermi_plus_log (data PbPb nominal; fitter lines 221-254, entry 459-465)
// Efficiency points mirror the data graph formation (RDFBasedHistFillingData.cxx:483-508
// + Utilities/HistFillUtils.h:49-82): per fine q·η range, ProjectionY("e") of the
// num/denom 2D over [bin_number(lo)+1, bin_number(hi)], then TGraphAsymmErrors::BayesDivide.
// (Caveat, same as data: BayesDivide on weighted hists uses effective entries.)
//
// Input : <sample dir>/mc_trig_eff_hists_<label>.root  (from FillMCTrigEffHists step 1)
// Output: <sample dir>/single_mu_effcy_pT_fit_mc.root
//   - TF1  f_mc_pt_vs_q_eta_<muplus|muminus>_<lo>_TO_<hi>   (pairToSuffix format)
//   - TGraphAsymmErrors g_mc_pt_vs_q_eta_<chg>_<lo>_TO_<hi> (the fitted points)
//   - TH2D h_mc_pt_vs_q_eta_ratio_<chg>                     (unfitted 2D fallback, num/denom)
// PNGs  : <sample dir>/mc_trig_eff_fit_plots/mc_trg_effcy_pT_fitting_<label>_<mu+|mu->.png
//         (one canvas per charge, all 10 fine q·η pads, log-x — data fitter layout)
//
// CORRECTED-MC study (corrected_mc = true): identical fit applied to the SF-corrected Step-1
// hists (mc_trig_eff_hists_<label><wp>_corrected.root) -> single_mu_effcy_pT_fit_mc_corrected<wp>.root.
// Object names inside the file are unchanged, so MCEffEvaluator loads ε_corr unmodified.
//
// Usage (from Analysis/RDFBasedHistFilling/):
//   root -b -l -q 'FitMCSinglesEffcy.cxx+("pp")'
//   root -b -l -q 'FitMCSinglesEffcy.cxx+("overlay")'
//   root -b -l -q 'FitMCSinglesEffcy.cxx+("pp_full", true, true)'   // corrected MC
// =============================================================================

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

using namespace std;

#include "../Utilities/bin_number.cxx"
#include "../Utilities/proj_range_to_suffix.cxx"
#include "CommonEffcyConfig.h"

namespace MCSinglesFit {

// same functional forms / init / limits / options as the data fitter (see header)
enum FittingMode { erf_plus_log, fermi_plus_log };

TF1* FitTurnOn(TGraphAsymmErrors* g, FittingMode mode, const std::string& fname, int& fit_status) {
    const double pT_min = 4;
    const double pT_max = 60;

    TF1* fTurnOn = nullptr;
    if (mode == erf_plus_log) {
        // data fitter fitTurnOnErfPlusLog (SingleMuEffcyPtTurnOnFitter.cxx:183-219)
        const std::string formula =
            "[2]*0.5*(1.0+TMath::Erf((x-[0])/(sqrt(2.0)*[1])))*(1.0+[3]*TMath::Log(1.0+(x-4.0)/4.0))";
        fTurnOn = new TF1(fname.c_str(), formula.c_str(), pT_min, pT_max);
        fTurnOn->SetParNames("mean", "sigma", "plateau", "corrCoef");
        fTurnOn->SetParameters(4.0, 2.0, 1.0, 0.);
        fTurnOn->SetParLimits(0, 0, 10);
        fTurnOn->SetParLimits(1, 0.01, 50);
        fTurnOn->SetParLimits(2, 0.5, 1.2);
        fTurnOn->SetParLimits(3, 0, 0.1);
    } else {
        // data fitter fitTurnOnFermiPlusLog (SingleMuEffcyPtTurnOnFitter.cxx:221-254)
        const std::string formula =
            "[0]/(1.0+TMath::Exp(([1]-x)/[2]))*(1.0+[3]*TMath::Log(1.0+(x-4.0)/4.0))";
        fTurnOn = new TF1(fname.c_str(), formula.c_str(), pT_min, pT_max);
        fTurnOn->SetParNames("normFermi", "pT0", "Delta", "corrCoef");
        fTurnOn->SetParameters(0.9, 4.0, 1.5, 0.);
        fTurnOn->SetParLimits(0, 0.6, 1.0);
        fTurnOn->SetParLimits(1, 2.5, 5.5);
        fTurnOn->SetParLimits(2, 0, 10);
        fTurnOn->SetParLimits(3, 0, 0.1);
    }

    // "QR" = the data fitter's options; the returned int is the fit status (0 = OK)
    fit_status = g->Fit(fTurnOn, "QR");
    return fTurnOn;
}

// ---------------------------------------------------------------------------------------------
// CORRECTED-MC efficiency graph (corrected_mc only).
//
// WHY NOT BayesDivide: the corrected numerator is the SF-re-weighted subset
// N = sum_fired w*SF, which EXCEEDS the denominator D = sum_all w in any bin where SF > 1.
// TGraphAsymmErrors::BayesDivide then rejects the pair ("passed TEfficiency objects do not have
// consistent bin contents") and returns an EMPTY graph -- and TGraph::Fit on an empty graph
// leaves the TF1 at its INITIAL parameters while still reporting status 0. That is a silent,
// badly wrong turn-on (seen at 3+ q.eta bins of the overlay before this was fixed), so the
// corrected graph is built explicitly instead.
//
// Central value: eff = N/D, exactly the estimator under study (eps_corr = eps_MC*<SF>).
// Error: the CONDITIONAL (binomial-correct) form, Var(N) = A - B with
//        A = sum_fired (w*SF)^2 and B = sum_fired (w*SF)^2 * eps_MC (booked by
//        FillMCTrigEffHists in corrected mode), e_eff = sqrt(A-B)/D.
//        BOTH binomial boundary cases (A-B <= 0 within round-off) get the "1/n rule" on the
//        effective denominator count n_eff = (D/e_D)^2:
//          k = n  (every effective entry fired) -> eff ~ 1, e ~ 1/n_eff;
//          k = 0  (none fired: A = B = 0)       -> eff = 0, e ~ 1/n_eff.
//        The k = 0 case is why the error is max(eff, 1)/n_eff and not eff/n_eff: an
//        eff = 0 +- 0 point is infinitely constraining and DESTROYS the fit -- it produced
//        chi2/ndf = 30627/30 and an unphysical eps(30 GeV) = 1.40 in pp_full q.eta [1.3,1.6) mu+
//        before this was fixed. (BayesDivide gives such a point a proper Bayesian upper error,
//        which is why the nominal path never saw this.)
// Zero-width bins (the data binning's duplicated 8.0 GeV edge) and empty denominators are skipped,
// exactly as BayesDivide drops them in the nominal path.
TGraphAsymmErrors* CorrectedEffGraph(const TH1* num, const TH1* den,
                                     const TH1* A, const TH1* B)
{
    auto* g = new TGraphAsymmErrors();
    int k = 0;
    for (int i = 1; i <= den->GetNbinsX(); ++i) {
        const double D = den->GetBinContent(i);
        const double bw = den->GetBinWidth(i);
        if (D <= 0. || bw <= 0.) continue;
        const double N = num->GetBinContent(i);
        const double a = A->GetBinContent(i), b = B->GetBinContent(i);
        double var = a - b;
        double e;
        if (var > 1e-6 * (a + b)) {
            e = std::sqrt(var) / D;
        } else {
            const double eD = den->GetBinError(i);
            const double neff = (eD > 0.) ? (D / eD) * (D / eD) : 1.0;
            e = (neff > 0.) ? std::max(N / D, 1.0) / neff : 0.;
        }
        g->SetPoint(k, den->GetBinCenter(i), N / D);
        g->SetPointError(k, 0.5 * bw, 0.5 * bw, e, e);
        ++k;
    }
    return g;
}

} // namespace MCSinglesFit

// =============================================================================
void FitMCSinglesEffcy(const std::string& sample = "pp", bool use_tight_wp = true,
                       bool corrected_mc = false, bool sf_closure = false) {
    using namespace MCSinglesFit;
    if (sf_closure && !corrected_mc) {
        std::cerr << "FitMCSinglesEffcy: sf_closure is a mode OF the corrected-MC path; "
                     "it needs corrected_mc = true" << std::endl;
        return;
    }

    std::string dir, label;
    FittingMode mode;
    if (sample == "pp") {
        dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/";
        label = "pp24";
        mode = erf_plus_log;    // data pp nominal
    } else if (sample == "pp_full") {
        // pp24 FULL sample: identical fit to "pp" (same pp turn-on shape); only the input dir
        // and the "pp24_full" label differ, so its outputs never clobber the TEST-sample fits.
        dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/";
        label = "pp24_full";
        mode = erf_plus_log;
    } else if (sample == "overlay") {
        dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/";
        label = "hijing_overlay_pbpb23";
        mode = fermi_plus_log;  // data PbPb nominal
    } else if (sample == "noovl") {
        // r17663 NO-OVERLAY diagnostic (R8, round 4): pp COLLISIONS -> the pp fit mode.
        // The comparison this sample exists for is against pp24 fullsim, so it must be
        // fitted with the same functional form; the reco conditions being PbPb-like does
        // not change the shape of a pp turn-on.
        dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_no_overlay_test_sample/";
        label = "r17663_no_overlay";
        mode = erf_plus_log;
    } else {
        std::cerr << "FitMCSinglesEffcy: sample must be \"pp\", \"pp_full\", \"overlay\" or "
                     "\"noovl\", got " << sample << std::endl;
        return;
    }

    // WP config (registry: Analysis/docs/muon_wp_registry.md): TIGHT nominal unsuffixed
    const std::string wp_suf = use_tight_wp ? "" : "_medium_wp";
    // CORRECTED-MC study (mc_trigger_efficiency.md round-7 contract item 5): fit the turn-ons of
    // the SF-corrected MC (numerator weighted by ε_data/ε_MC) with the SAME functional form,
    // parameter init/limits and fit options as the nominal MC, so the corrected and the nominal
    // fits differ ONLY by the corrected numerator. The OBJECT names inside the file are
    // deliberately UNCHANGED (f_mc_pt_vs_q_eta_*, g_mc_pt_vs_q_eta_*, h_mc_pt_vs_q_eta_ratio_*)
    // so that MCEffEvaluator can load ε_corr with no special casing; only the FILE name differs,
    // which is what keeps the nominal fit file from ever being overwritten.
    // sf_closure: the SF ≡ 1 closure variant of the corrected chain (see FillMCTrigEffHists).
    const std::string corr_suf = corrected_mc ? (sf_closure ? "_corrected_sfclosure" : "_corrected")
                                              : "";
    const std::string infile_name  = dir + "mc_trig_eff_hists_" + label + wp_suf + corr_suf + ".root";
    const std::string outfile_name = dir + "single_mu_effcy_pT_fit_mc" + corr_suf + wp_suf + ".root";
    const std::string plot_dir = dir + "mc_trig_eff_fit_plots/";
    gSystem->mkdir(plot_dir.c_str(), kTRUE);

    std::cout << "FitMCSinglesEffcy: sample=" << sample << " (" << label << "), mode="
              << (mode == erf_plus_log ? "erf_plus_log" : "fermi_plus_log")
              << ", corrected_mc=" << corrected_mc << std::endl;
    std::cout << "  in : " << infile_name << "\n  out: " << outfile_name << std::endl;

    TFile* fin = TFile::Open(infile_name.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "FitMCSinglesEffcy: cannot open " << infile_name << std::endl;
        return;
    }
    TFile* fout = TFile::Open(outfile_name.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "FitMCSinglesEffcy: cannot open output " << outfile_name << std::endl;
        return;
    }

    static const CommonEffcyConfig cfg{};
    gStyle->SetOptStat(0);

    int n_fits = 0, n_failed = 0;

    for (const std::string chg : {"muplus", "muminus"}) {
        TH2D* h_num = dynamic_cast<TH2D*>(fin->Get(("h_mc_pt_vs_q_eta_num_" + chg).c_str()));
        TH2D* h_den = dynamic_cast<TH2D*>(fin->Get(("h_mc_pt_vs_q_eta_denom_" + chg).c_str()));
        if (!h_num || !h_den) {
            std::cerr << "FitMCSinglesEffcy: missing 2D hists for " << chg << " in " << infile_name << std::endl;
            continue;
        }
        // corrected mode: the conditional-error terms of the corrected numerator (see
        // CorrectedEffGraph). Their absence is FATAL, not a reason to fall back to BayesDivide --
        // that path is undefined for a re-weighted numerator and fails silently.
        TH2D* h_A = nullptr;
        TH2D* h_B = nullptr;
        if (corrected_mc) {
            h_A = dynamic_cast<TH2D*>(fin->Get(("h_mc_pt_vs_q_eta_errA_" + chg).c_str()));
            h_B = dynamic_cast<TH2D*>(fin->Get(("h_mc_pt_vs_q_eta_errB_" + chg).c_str()));
            if (!h_A || !h_B) {
                std::cerr << "FitMCSinglesEffcy: corrected mode needs h_mc_pt_vs_q_eta_errA/errB_"
                          << chg << " in " << infile_name
                          << " -- rerun FillMCTrigEffHists with corrected_mc = true" << std::endl;
                fout->Close(); fin->Close();
                return;
            }
        }

        // unfitted 2D ratio fallback (gap q·η regions in Step 3)
        fout->cd();
        TH2D* h_ratio = static_cast<TH2D*>(h_num->Clone(("h_mc_pt_vs_q_eta_ratio_" + chg).c_str()));
        h_ratio->Divide(h_den);
        h_ratio->Write();

        // q·η axis edges from the histogram itself (mirrors RDFBasedHistFillingData.cxx:483-489)
        std::vector<double> q_eta_edges;
        const TAxis* xaxis = h_den->GetXaxis();
        const double* arr = xaxis->GetXbins()->GetArray();
        q_eta_edges.assign(arr, arr + xaxis->GetNbins() + 1);

        // data-fitter canvas layout: 4x3, log-x pads (SingleMuEffcyPtTurnOnFitter.cxx:310-317)
        const std::string chg_label = (chg == "muplus") ? "mu+" : "mu-";
        TCanvas* c = new TCanvas(("c_mc_" + chg).c_str(), "MC Trigger Turn-on Curves", 1500, 1000);
        c->Divide(4, 3);

        int idx = 0;
        for (const auto& range : cfg.q_eta_proj_ranges_fine_excl_gap) {
            const std::string q_eta_suffix = pairToSuffix(range);

            // projection bin range, exactly as the data graph maker (Data.cxx:491-493)
            const int bin_first = bin_number(range.first, q_eta_edges) + 1;
            const int bin_last  = bin_number(range.second, q_eta_edges);

            std::unique_ptr<TH1D> h_num1D(h_num->ProjectionY(
                Form("%s_py_%s", h_num->GetName(), q_eta_suffix.c_str()), bin_first, bin_last, "e"));
            std::unique_ptr<TH1D> h_den1D(h_den->ProjectionY(
                Form("%s_py_%s", h_den->GetName(), q_eta_suffix.c_str()), bin_first, bin_last, "e"));

            // Graph builder. NOMINAL -- and the SF ≡ 1 CLOSURE, whose numerator is bit-identical
            // to the nominal one so num <= denom is guaranteed -- use BayesDivide, so that the
            // closure exercises exactly the nominal estimator and its Step-3/4 output must
            // reproduce the nominal output bin-by-bin. Only the genuine corrected numerator needs
            // CorrectedEffGraph (see its header).
            TGraphAsymmErrors* g = nullptr;
            if (!corrected_mc || sf_closure) {
                g = new TGraphAsymmErrors();
                g->BayesDivide(h_num1D.get(), h_den1D.get());
            } else {
                std::unique_ptr<TH1D> h_A1D(h_A->ProjectionY(
                    Form("%s_py_%s", h_A->GetName(), q_eta_suffix.c_str()), bin_first, bin_last, "e"));
                std::unique_ptr<TH1D> h_B1D(h_B->ProjectionY(
                    Form("%s_py_%s", h_B->GetName(), q_eta_suffix.c_str()), bin_first, bin_last, "e"));
                g = CorrectedEffGraph(h_num1D.get(), h_den1D.get(), h_A1D.get(), h_B1D.get());
            }
            g->SetName(("g_mc_pt_vs_q_eta_" + chg + "_" + q_eta_suffix).c_str());

            int fit_status = -1;
            TF1* fit = FitTurnOn(g, mode, "f_mc_pt_vs_q_eta_" + chg + "_" + q_eta_suffix, fit_status);
            ++n_fits;
            const double plateau_val = fit->Eval(30.0);
            std::cout << "  " << fit->GetName() << ": status=" << fit_status
                      << ", Npts=" << g->GetN()
                      << ", chi2/ndf=" << fit->GetChisquare() << "/" << fit->GetNDF()
                      << ", eps(30 GeV)=" << plateau_val << std::endl;
            // An EMPTY graph makes TGraph::Fit a no-op that still returns status 0 and leaves the
            // TF1 at its initial parameters -- exactly the silent failure this fitter must never
            // ship. Treat "no points" and "no degrees of freedom" as fit failures.
            if (fit_status != 0 || g->GetN() == 0 || fit->GetNDF() <= 0) {
                ++n_failed;
                std::cerr << "  WARNING: fit FAILED (status " << fit_status << ", Npts="
                          << g->GetN() << ", ndf=" << fit->GetNDF() << "): "
                          << fit->GetName() << std::endl;
            }

            fout->cd();
            fit->Write();
            g->Write();

            // overlay pad (data-fitter styling)
            c->cd(idx + 1);
            gPad->SetLogx();
            g->SetMarkerColor(kBlack);
            g->SetLineColor(kBlack);
            g->SetMarkerStyle(20);
            g->SetMarkerSize(0.9);
            fit->SetLineColor(kRed);
            fit->SetLineWidth(2);
            g->GetXaxis()->SetTitle("p_{T} [GeV]");
            g->GetYaxis()->SetTitle("#epsilon");
            g->GetYaxis()->SetRangeUser(0, 1.1);
            g->GetXaxis()->SetLimits(4.0, 60.0);
            g->Draw("AP");
            fit->Draw("SAME");

            TLegend* leg = new TLegend(0.35, 0.25, 0.88, 0.5);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(g, ("MC mu4, " + std::string(chg == "muplus" ? "#mu^{+}" : "#mu^{-}")).c_str(), "lp");
            leg->AddEntry(static_cast<TObject*>(nullptr), pairToLegendLabel(range).c_str(), "");
            leg->Draw("SAME");

            ++idx;
        }

        // NOMINAL png name left byte-identical (no WP token -- a pre-existing quirk: the Medium
        // run overwrites the Tight png). The CORRECTED pngs carry BOTH tokens so that they
        // neither clobber the nominal ones nor each other across working points.
        const std::string png_tag = corrected_mc ? (corr_suf + wp_suf) : std::string("");
        c->SaveAs(Form("%smc_trg_effcy_pT_fitting_%s%s_%s.png", plot_dir.c_str(), label.c_str(),
                       png_tag.c_str(), chg_label.c_str()));
    }

    fout->Close();
    fin->Close();
    std::cout << "FitMCSinglesEffcy: " << n_fits << " fits, " << n_failed << " failed. Output: "
              << outfile_name << " + PNGs in " << plot_dir << std::endl;
}
