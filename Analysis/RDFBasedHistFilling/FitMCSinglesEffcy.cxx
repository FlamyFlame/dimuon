// =============================================================================
// FitMCSinglesEffcy.cxx
//
// MC single-muon mu4 turn-on fits (docs/tracking/mc_trigger_efficiency.md §3.1),
// consumed by FillMCTrigEffHists Step 3 (§3.3 inverse weights).
//
// Mirrors the DATA fitter Analysis/SingleMuEffcyPtTurnOnFitter.cxx in the fit
// ingredients (same functional forms, parameter init/limits, "QR" fit option, and the
// same TFormula-string TF1 so the analytic formula persists on Write), with ONE
// DELIBERATE DIFFERENCE:
//
//   FIT RANGE. MC is [4.5,60] with the log/linear pivot at 4.5; DATA is [4,60] with the
//   pivot at 4.0. This is NOT drift -- it is the two samples' populations differing on
//   purpose (docs/tracking/mu_pt45_gap125_pairpt9_adoption.md D8/D9/D11). MC has no muons
//   below 4.5 GeV at all: the MC trigger efficiency is truth-seeded and reco-matched, not
//   tag-and-probe, so there is no probe population to widen, and its deliverable eps_dR is
//   pair-level and MUST be measured on the 4.5 GeV analysis population. The DATA fitter
//   keeps 4.0 because its tag-and-probe eps^nc is a per-muon efficiency whose measurement
//   population may be wider than the one it is applied to.
//   Both are evaluated only above 4.5 GeV, where the analysis lives, and pT_bins_8 forces
//   4.5 to be a bin edge so the two maps stay bit-identical there.
//
// Everything else is a true mirror:
//   pp      -> erf_plus_log  (data pp nominal;   fitter lines 183-219, entry 451-457)
//   overlay -> fermi_plus_log (data PbPb nominal; fitter lines 221-254, entry 459-465)
// Efficiency points mirror the data graph formation (RDFBasedHistFillingData.cxx:483-508
// + Utilities/HistFillUtils.h:49-82): per fine q·η range, ProjectionY("e") of the
// num/denom 2D over [bin_number(lo)+1, bin_number(hi)], then TGraphAsymmErrors::BayesDivide.
// (Caveat, same as data: BayesDivide on weighted hists uses effective entries.)
//
// Input : <sample dir>/mc_trig_eff/hists/mc_trig_eff_hists_<label>.root  (FillMCTrigEffHists step 1)
// Output: <sample dir>/mc_trig_eff/singles_fits/single_mu_effcy_pT_fit_mc_<label>.root
//         (FullSimMCSinglesFitFile; the label is part of the basename since 2026-09-16 --
//          mc_trigger_efficiency.md RW 3c -- so a wrong directory can no longer overwrite a
//          sibling sample's fit silently)
//   - TF1  f_mc_pt_vs_q_eta_<muplus|muminus>_<lo>_TO_<hi>   (pairToSuffix format)
//   - TGraphAsymmErrors g_mc_pt_vs_q_eta_<chg>_<lo>_TO_<hi> (the fitted points)
//   - TH2D h_mc_pt_vs_q_eta_ratio_<chg>                     (unfitted 2D fallback, num/denom)
// No PNGs (retired 2026-09-16, user decision): the fit is drawn ONLY by the Step-1 panel
// plot_mc_trig_eff.cxx -> step1_singles_data_mc/step1_eff_pt_in_q_eta_bins_<charge>.png, which
// overlays the same points and the same TF1 on the data reference. The former
// mc_trig_eff_fit_plots/ canvases were a strict subset of it.
//
// CORRECTED-MC study (corrected_mc = true): identical fit applied to the SF-corrected Step-1
// hists (mc_trig_eff_hists_<label><wp>_corrected.root) -> single_mu_effcy_pT_fit_mc_<label>_corrected<wp>.root.
// Object names inside the file are unchanged, so MCEffEvaluator loads ε_corr unmodified.
//
// overlay_year (LAST argument, default 24): which HIJING-overlay TEST production "overlay" means
// (24 = Pb+Pb 2024 conditions, 23 = ..._pbpb23/); see dr_correction_sample_cfg.h.
//
// Usage (from Analysis/RDFBasedHistFilling/):
//   root -b -l -q 'FitMCSinglesEffcy.cxx+("pp")'
//   root -b -l -q 'FitMCSinglesEffcy.cxx+("overlay")'
//   root -b -l -q 'FitMCSinglesEffcy.cxx+("overlay", true, false, false, 23)'   // pbpb23 overlay
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

#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TString.h>   // Form(), used to compose the log pivot from pT_min
#include <TSystem.h>

using namespace std;

#include "../Utilities/bin_number.cxx"
#include "../Utilities/proj_range_to_suffix.cxx"
#include "CommonEffcyConfig.h"
#include "../plotting_codes/trig_effcy/mc_based/dr_correction_sample_cfg.h"   // sample table + product files

namespace MCSinglesFit {

// same functional forms / init / limits / options as the data fitter (see header)
enum FittingMode { erf_plus_log, fermi_plus_log };

TF1* FitTurnOn(TGraphAsymmErrors* g, FittingMode mode, const std::string& fname, int& fit_status) {
    // Low edge 4.5 GeV, per D9: MC has no muons below 4.5 at all. This does NOT mirror the data
    // fitter, which deliberately stays at 4.0 -- see this file's header for the ONE DELIBERATE
    // DIFFERENCE (D11) and why it is not drift.
    // The pivot of the log correction term IS pT_min and is COMPOSED from it
    // rather than retyped -- the two used to be independent "4.0" literals in both fitters,
    // which is precisely how a range change silently leaves the pivot behind.
    const double pT_min = 4.5;
    const double pT_max = 60;
    const std::string pivot = Form("%g", pT_min);

    // Parameter inits/limits anchored near 4 GeV (erf mean, Fermi pT0 and its [2.5,5.5]
    // limits) are PHYSICS PRIORS on where the mu4 turn-on sits -- a trigger/detector
    // property, unchanged by where the offline analysis cuts -- so they are deliberately NOT
    // moved with the fit range, in both fitters.
    TF1* fTurnOn = nullptr;
    if (mode == erf_plus_log) {
        // data fitter fitTurnOnErfPlusLog (SingleMuEffcyPtTurnOnFitter.cxx:193-238)
        const std::string formula =
            "[2]*0.5*(1.0+TMath::Erf((x-[0])/(sqrt(2.0)*[1])))*(1.0+[3]*TMath::Log(1.0+(x-"
            + pivot + ")/" + pivot + "))";
        fTurnOn = new TF1(fname.c_str(), formula.c_str(), pT_min, pT_max);
        fTurnOn->SetParNames("mean", "sigma", "plateau", "corrCoef");
        fTurnOn->SetParameters(4.0, 2.0, 1.0, 0.);
        fTurnOn->SetParLimits(0, 0, 10);
        fTurnOn->SetParLimits(1, 0.01, 50);
        fTurnOn->SetParLimits(2, 0.5, 1.2);
        fTurnOn->SetParLimits(3, 0, 0.1);
    } else {
        // data fitter fitTurnOnFermiPlusLog (SingleMuEffcyPtTurnOnFitter.cxx:240-280)
        const std::string formula =
            "[0]/(1.0+TMath::Exp(([1]-x)/[2]))*(1.0+[3]*TMath::Log(1.0+(x-"
            + pivot + ")/" + pivot + "))";
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
                       bool corrected_mc = false, bool sf_closure = false,
                       int overlay_year = 24) {
    using namespace MCSinglesFit;
    if (sf_closure && !corrected_mc) {
        std::cerr << "FitMCSinglesEffcy: sf_closure is a mode OF the corrected-MC path; "
                     "it needs corrected_mc = true" << std::endl;
        return;
    }

    // Sample identity (directory, label) from the one table the whole chain shares. Throws on an
    // unknown key.
    const DrCorrSample cfg = GetDrCorrSample(sample, use_tight_wp, overlay_year);
    const std::string label = cfg.mc_label;
    // Fit form = the data-side nominal of the collision system the sample simulates:
    //   pp collisions (pp, pp_full, and the r17663 no-overlay diagnostic -- pp COLLISIONS
    //   reconstructed with PbPb conditions; its whole purpose is the comparison against pp24
    //   fullsim, so it must be fitted with the same form) -> erf_plus_log (data pp nominal);
    //   HIJING overlay (Pb+Pb)                             -> fermi_plus_log (data PbPb nominal).
    const FittingMode mode = (sample == "overlay") ? fermi_plus_log : erf_plus_log;
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
    const std::string infile_name  = DrCorrHistFile(cfg, use_tight_wp, corr_suf);
    const std::string outfile_name = DrCorrSinglesFitFile(cfg, use_tight_wp, corr_suf);
    gSystem->mkdir(gSystem->DirName(outfile_name.c_str()), kTRUE);

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

    static const CommonEffcyConfig eff_cfg{};

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

        for (const auto& range : eff_cfg.q_eta_proj_ranges_coarse_incl_gap) {
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
        }
    }

    fout->Close();
    fin->Close();
    std::cout << "FitMCSinglesEffcy: " << n_fits << " fits, " << n_failed << " failed. Output: "
              << outfile_name << std::endl;
}
