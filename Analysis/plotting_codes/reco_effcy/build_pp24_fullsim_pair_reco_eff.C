// =================================================================================================
// build_pp24_fullsim_pair_reco_eff.C
//
// Persist the single-b PAIR reconstruction efficiency eps_reco(pair pT, pair eta, dR) that the
// pp24 cross-section applies, from the pp24-condition Pythia fullsim RDF hist-filling output.
//
// WHY THIS FILE EXISTS. Until now the fullsim reco-efficiency chain produced PNGs only -- there was
// no persisted efficiency product anywhere, so the cross-section had to fall back on the Run-2
// single-muon PLACEHOLDER (eps_1*eps_2, no dR dependence;
// docs/tracking/reco_eff_placeholder_run2.md). This macro writes the ROOT product the
// cross-section reads through Utilities/PairRecoEffEvaluator.h.
//
// DEFINITION (docs/tracking/pp24_crossx_rerun_2026_08.md Physics Procedure 3d):
//
//   eps_reco(cell) = N[ single-b OS pair, both muons truth-matched, pair passes the WP, and the
//                       RECO pair passes the signal region incl. the fiducial gap cut on RECO q*eta ]
//                  / N[ single-b OS pair whose TRUTH pair passes the signal region
//                       incl. the fiducial gap cut on TRUTH q*eta ]
//
// Numerator and denominator are the `_single_b_pass_{tight,medium}_and_signal_truth_and_reco` and
// `_single_b_pass_signal_truth` filters of RDFBasedHistFillingPythiaFullsim, binned in TRUTH
// kinematics on the CANONICAL coarse axes and weighted by the MC event weight. Because the gap cut
// sits in BOTH legs, eps_reco is a FIDUCIAL efficiency: it does NOT contain the truth-level gap
// acceptance eps_acc (muon_gap_cuts_acceptance.md F12), which stays a separate factor.
//
// TWO FALLBACK LEVELS ARE WRITTEN ALONGSIDE THE 3D MAP, and they are not optional. The three
// variables are strongly correlated inside the signal region (dR ~< 2 m_uu / pT^pair, so high
// pair pT forces small dR): about half of the 3D cells have NO denominator at all (measured: 143 of 288 empty, plus one
// cell with a denominator but no reconstructed pair, so 144 of 288 deliver no efficiency). A pair landing
// in one of those must not be left uncorrected -- w_reco = 1 would be a silent, one-sided bias.
// So the evaluator falls back, per pair and counted:
//     3D cell  ->  the dR-INTEGRATED efficiency of the same (pair pT, pair eta) cell
//              ->  the INCLUSIVE efficiency of the whole signal region.
//
// Usage:  root -l -b -q 'build_pp24_fullsim_pair_reco_eff.C+()'
//         root -l -b -q 'build_pp24_fullsim_pair_reco_eff.C+(false)'   // TEST sample instead
// =================================================================================================

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNamed.h>
#include <TParameter.h>

namespace {

const char* kBase =
    "h_truth_dr_effcy_vs_truth_pair_eta_coarse_vs_truth_pair_pt_coarse_single_b";

TH3D* Get3D(TFile* f, const std::string& name)
{
    auto* h = dynamic_cast<TH3D*>(f->Get(name.c_str()));
    if (!h) {
        std::cerr << "build_pp24_fullsim_pair_reco_eff: MISSING " << name << " in "
                  << f->GetName() << "\n  -- rerun the fullsim RDF hist filling "
                     "(pipelines/pipeline_pythia_fullsim_pp.sh, stage 5) after the "
                     "2026-08-17 binning addition." << std::endl;
        return nullptr;
    }
    return h;
}

}  // namespace

int build_pp24_fullsim_pair_reco_eff(bool use_full_sample = true)
{
    const std::string dir = use_full_sample
        ? "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/"
        : "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/";
    const std::string sfx = use_full_sample ? "_full" : "";
    const std::string in_path  = dir + "histograms_pythia_fullsim_pp24_no_data_resonance_cuts"
                               + sfx + ".root";
    const std::string out_path = dir + "pair_reco_eff_pp24" + sfx + ".root";

    std::unique_ptr<TFile> fin(TFile::Open(in_path.c_str(), "READ"));
    if (!fin || fin->IsZombie()) {
        std::cerr << "build_pp24_fullsim_pair_reco_eff: cannot open " << in_path << std::endl;
        return 1;
    }

    TH3D* h_den = Get3D(fin.get(), std::string(kBase) + "_pass_signal_truth");
    if (!h_den) return 2;

    std::unique_ptr<TFile> fout(TFile::Open(out_path.c_str(), "RECREATE"));

    for (const std::string wp : {std::string("tight"), std::string("medium")}) {
        TH3D* h_num = Get3D(fin.get(),
                            std::string(kBase) + "_pass_" + wp + "_and_signal_truth_and_reco");
        if (!h_num) { fout->Close(); return 3; }

        // ---- 3D map ----------------------------------------------------------------------
        auto* eff = static_cast<TH3D*>(h_num->Clone(("h_pair_reco_eff_" + wp).c_str()));
        eff->SetDirectory(nullptr);
        eff->Divide(h_num, h_den, 1.0, 1.0, "B");
        eff->SetTitle(("single-b pair reco efficiency, " + wp
                     + " WP;p_{T}^{pair,truth} [GeV];#eta^{pair,truth};#DeltaR^{truth}").c_str());

        // ---- dR-integrated 2D fallback ---------------------------------------------------
        h_num->SetName(("num3d_" + wp).c_str());
        h_den->SetName("den3d");
        auto* n2 = static_cast<TH2D*>(h_num->Project3D("yx"));   // -> "num3d_<wp>_yx"
        auto* d2 = static_cast<TH2D*>(h_den->Project3D("yx"));   // -> "den3d_yx"
        auto* eff2 = static_cast<TH2D*>(n2->Clone(("h_pair_reco_eff_" + wp + "_dr_integrated").c_str()));
        eff2->SetDirectory(nullptr);
        eff2->Divide(n2, d2, 1.0, 1.0, "B");
        eff2->SetTitle(("single-b pair reco efficiency, dR-integrated, " + wp
                      + " WP;p_{T}^{pair,truth} [GeV];#eta^{pair,truth}").c_str());

        // ---- inclusive scalar fallback ---------------------------------------------------
        const double N = h_num->Integral(), D = h_den->Integral();
        const double eff_incl = (D > 0.) ? N / D : -1.0;

        // ---- census, so the placeholder fraction is never invisible ----------------------
        int n_cells = 0, n_empty = 0, n_empty2 = 0, n_cells2 = 0;
        double w_in_empty = 0., w_tot = 0.;
        for (int ix = 1; ix <= eff->GetNbinsX(); ++ix)
            for (int iy = 1; iy <= eff->GetNbinsY(); ++iy) {
                for (int iz = 1; iz <= eff->GetNbinsZ(); ++iz) {
                    ++n_cells;
                    const double d = h_den->GetBinContent(ix, iy, iz);
                    w_tot += d;
                    // "No measure" = the cell delivers no usable efficiency: either the truth
                    // denominator is empty, or it is populated but nothing was reconstructed
                    // there. Both route to the dR-integrated fallback; the SECOND kind is the
                    // one that costs signal, so its denominator weight is tracked separately.
                    if (d <= 0.) { ++n_empty; }
                    else if (eff->GetBinContent(ix, iy, iz) <= 0.) { ++n_empty; w_in_empty += d; }
                }
                ++n_cells2;
                if (d2->GetBinContent(ix, iy) <= 0.) ++n_empty2;
            }

        std::cout << "\n=== eps_reco, " << wp << " WP ===\n"
                  << "  3D cells                 : " << n_cells << " (" << eff->GetNbinsX()
                  << " pair pT x " << eff->GetNbinsY() << " pair eta x " << eff->GetNbinsZ()
                  << " dR)\n"
                  << "  3D cells with no measure : " << n_empty << "  ("
                  << (n_cells ? 100.0 * n_empty / n_cells : 0.0) << " % of cells) -> dR-integrated"
                     " fallback\n"
                  << "  2D (pT,eta) cells empty  : " << n_empty2 << " of " << n_cells2
                  << " -> inclusive fallback\n"
                  << "  MC denominator weight in an empty 3D cell : "
                  << (w_tot > 0. ? 100.0 * w_in_empty / w_tot : 0.0) << " %\n"
                  << "  INCLUSIVE eps_reco       : " << eff_incl << std::endl;

        fout->cd();
        eff->Write();
        eff2->Write();
        TParameter<double> p_incl(("pair_reco_eff_inclusive_" + wp).c_str(), eff_incl);
        p_incl.Write();
        // Keep the raw ingredients in the same file: a consumer that wants a different error
        // treatment must not have to guess which histograms these came from.
        h_num->Clone(("h_pair_reco_eff_num_" + wp).c_str())->Write();
        n2->Clone(("h_pair_reco_eff_num2d_" + wp).c_str())->Write();
    }

    fout->cd();
    h_den->Clone("h_pair_reco_eff_denom")->Write();
    static_cast<TH2D*>(h_den->Project3D("yx"))->Clone("h_pair_reco_eff_denom2d")->Write();
    TNamed prov("provenance",
                ("single-b OS pair reco efficiency; fiducial (gap cut on BOTH the truth and the "
                 "reco leg, ParamsSet::single_mu_fiducial_gap_cuts); source " + in_path
                 + "; axes = ParamsSet::pair_pt_coarse_bins x "
                   "CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap x "
                   "RDFBasedHistFillingPythia::dr_bins_edges_for_reco_effcy; does NOT contain "
                   "eps_acc").c_str());
    prov.Write();
    fout->Close();

    std::cout << "\nWrote " << out_path << std::endl;
    return 0;
}
