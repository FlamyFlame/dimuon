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
// Usage (from Analysis/RDFBasedHistFilling/):
//   root -b -l -q 'FitMCSinglesEffcy.cxx+("pp")'
//   root -b -l -q 'FitMCSinglesEffcy.cxx+("overlay")'
// =============================================================================

#include <iostream>
#include <map>
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

} // namespace MCSinglesFit

// =============================================================================
void FitMCSinglesEffcy(const std::string& sample = "pp", bool use_tight_wp = true) {
    using namespace MCSinglesFit;

    std::string dir, label;
    FittingMode mode;
    if (sample == "pp") {
        dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/";
        label = "pp24";
        mode = erf_plus_log;    // data pp nominal
    } else if (sample == "overlay") {
        dir = "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/";
        label = "hijing_overlay_pbpb23";
        mode = fermi_plus_log;  // data PbPb nominal
    } else {
        std::cerr << "FitMCSinglesEffcy: sample must be \"pp\" or \"overlay\", got " << sample << std::endl;
        return;
    }

    // WP config (registry: Analysis/docs/muon_wp_registry.md): TIGHT nominal unsuffixed
    const std::string wp_suf = use_tight_wp ? "" : "_medium_wp";
    const std::string infile_name = dir + "mc_trig_eff_hists_" + label + wp_suf + ".root";
    const std::string outfile_name = dir + "single_mu_effcy_pT_fit_mc" + wp_suf + ".root";
    const std::string plot_dir = dir + "mc_trig_eff_fit_plots/";
    gSystem->mkdir(plot_dir.c_str(), kTRUE);

    std::cout << "FitMCSinglesEffcy: sample=" << sample << " (" << label << "), mode="
              << (mode == erf_plus_log ? "erf_plus_log" : "fermi_plus_log") << std::endl;

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

            auto* g = new TGraphAsymmErrors();
            g->BayesDivide(h_num1D.get(), h_den1D.get());
            g->SetName(("g_mc_pt_vs_q_eta_" + chg + "_" + q_eta_suffix).c_str());

            int fit_status = -1;
            TF1* fit = FitTurnOn(g, mode, "f_mc_pt_vs_q_eta_" + chg + "_" + q_eta_suffix, fit_status);
            ++n_fits;
            const double plateau_val = fit->Eval(30.0);
            std::cout << "  " << fit->GetName() << ": status=" << fit_status
                      << ", chi2/ndf=" << fit->GetChisquare() << "/" << fit->GetNDF()
                      << ", eps(30 GeV)=" << plateau_val << std::endl;
            if (fit_status != 0) {
                ++n_failed;
                std::cerr << "  WARNING: fit FAILED (status " << fit_status << "): " << fit->GetName() << std::endl;
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

        c->SaveAs(Form("%smc_trg_effcy_pT_fitting_%s_%s.png", plot_dir.c_str(), label.c_str(), chg_label.c_str()));
    }

    fout->Close();
    fin->Close();
    std::cout << "FitMCSinglesEffcy: " << n_fits << " fits, " << n_failed << " failed. Output: "
              << outfile_name << " + PNGs in " << plot_dir << std::endl;
}
