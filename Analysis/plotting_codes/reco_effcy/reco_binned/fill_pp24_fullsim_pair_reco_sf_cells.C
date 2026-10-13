// =================================================================================================
// fill_pp24_fullsim_pair_reco_sf_cells.C
//
// SF-weighted (pair pT, eta^pair, m_mumu) cell ingredients of the data-corrected reco-binned pair
// reco efficiency, and the NO-GAP reco leg used to scale the Pb+Pb Run-2 placeholder
// (docs/tracking/crossx_raa_reco_binned_propagation_20261013.md, Physics Procedure §2, §3a):
//
//   eps^{reco,SF}_j = sum_{R_j} w * SF_1 * SF_2 / N[T_j]       SF_i = m_i.eff_sf_<wp>
//   eps^{nogap}_p   = sum_{R^nogap_p} w / N[T_p]               (no single-muon q*eta gap windows)
//
// Input : the ntuple-processing output muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_full.root,
//         tree muon_pair_tree_sign2. It is the SAME truth-seeded pair sample as the nominal `_full` tree
//         (PythiaFullSimExtras.c ProcessEventFullsim: store_mc_trigger only adds trigger fields and the
//         per-muon SFs), the only one that carries eff_sf_<wp>. Set identity is ENFORCED: the unweighted
//         N[T_j], N[R_j], N[T_j and R_j] filled here must equal the migration file
//         (fill_pp24_fullsim_pair_bin_migration.C, built from `_full`) cell by cell, else this THROWS.
// Cuts  : identical expressions to fill_pp24_fullsim_pair_bin_migration.C (RDF class builders).
//
// Output: reco_eff/pair_reco_sf_cells_pp24_full.root
//   h_cell_minv_T, h_cell_minv_{R,TR}_<wp> (TH1D), h_cell_minv_M_<wp> (TH2D)        unweighted (identity check)
//   h_cell_minv_{R,TR}_sf_<wp>, h_cell_minv_M_sf_<wp>                              weight w*SF_1*SF_2 on R
//   h_pair_pt_T (truth pT, T), h_pair_pt_Rnogap_<wp> (reco pT, R without the q*eta gap windows)
//   (h_cell_minv_T is also stored as h_cell_minv_T so reco_binned_cells.h can read "sf_<wp>" sets.)
//
// Usage (from Analysis/plotting_codes/reco_effcy/reco_binned/):
//   root -l -b -q 'fill_pp24_fullsim_pair_reco_sf_cells.C+'
// =================================================================================================

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <TAxis.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TSystem.h>

#include "../../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../../RDFBasedHistFilling/RDFBasedHistFillingPythia.h"
#include "../../../Utilities/PairRecoEffRecoBinnedEvaluator.h"

namespace {

struct BinIdx {
    TAxis ax;
    explicit BinIdx(const std::vector<double>& e) : ax((int)e.size() - 1, e.data()) {}
    int operator()(double v) const
    {
        const int b = ax.FindFixBin(v);
        return (b >= 1 && b <= ax.GetNbins()) ? b - 1 : -1;
    }
};

void CompareIdentical(const TH1* mine, TFile* fmig, const std::string& name)
{
    auto* ref = dynamic_cast<TH1*>(fmig->Get(name.c_str()));
    if (!ref) throw std::runtime_error("identity check: migration file has no " + name);
    if (ref->GetNcells() != mine->GetNcells()) throw std::runtime_error("identity check: cell count differs for " + name);
    double max_rel = 0;
    for (int b = 0; b < ref->GetNcells(); ++b) {
        const double a = ref->GetBinContent(b), c = mine->GetBinContent(b);
        if (a == 0 && c == 0) continue;
        max_rel = std::max(max_rel, std::fabs(a - c) / std::max(std::fabs(a), std::fabs(c)));
    }
    std::cout << "[identity] " << name << " (mc_trig tree) vs migration file (nominal tree): max |rel diff| = "
              << max_rel << "\n";
    if (max_rel > 1e-6)
        throw std::runtime_error("identity check FAILED for " + name
                                 + ": the mc_trig pair tree is NOT the nominal pair sample -- do not use its SFs");
}

}  // namespace

void fill_pp24_fullsim_pair_reco_sf_cells()
{
    using EV = PairRecoEffRecoBinnedEvaluator;
    const std::vector<double> ept = EV::PtEdges(), eeta = EV::EtaEdges(), em = EV::MinvEdges();
    const int npt = ept.size() - 1, neta = eeta.size() - 1, nm = em.size() - 1;
    const int ncell = npt * neta * nm;

    const std::string sample_dir = FullSimSampleInputDir(FullSimSampleType::pp, /*is_test_sample=*/false);
    const std::string in_file  = sample_dir + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_full.root";
    const std::string reco_dir = FullSimRecoEffDir(sample_dir);
    const std::string mig_file = reco_dir + "pair_bin_migration_pp24_full.root";
    const std::string out_file = reco_dir + "pair_reco_sf_cells_pp24_full.root";

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df0("muon_pair_tree_sign2", in_file);

    const BinIdx bpt(ept), beta(eeta), bm(em);
    // Same strings as fill_pp24_fullsim_pair_bin_migration.C.
    const std::string sig_truth = ParamsSet::SignalMinvCutExpr("truth_minv") + " && "
                                + ParamsSet::SignalPairPtCutExpr("truth_pair_pt");
    const std::string sig_reco  = "(m1.reco_match && m2.reco_match) ? (" + ParamsSet::SignalMinvCutExpr("minv")
                                + " && " + ParamsSet::SignalPairPtCutExpr("pair_pt") + " && "
                                + RDFBasedHistFillingPythiaFullsim::RecoGapCutExpr() + ") : false";
    // NO-GAP reco leg: the same minus the single-muon q*eta windows; |eta^pair| < 2.2 kept (= the axis).
    const std::string sig_reco_nogap = "(m1.reco_match && m2.reco_match) ? (" + ParamsSet::SignalMinvCutExpr("minv")
                                + " && " + ParamsSet::SignalPairPtCutExpr("pair_pt") + " && "
                                + ParamsSet::PairFiducialEtaCutExpr("pair_eta") + ") : false";

    auto df = df0
        .Define("pass_signal_truth", sig_truth)
        .Define("pass_signal_reco",  sig_reco)
        .Define("pass_signal_reco_nogap", sig_reco_nogap)
        .Define("t_pt",  [bpt](float v)  { return bpt(v);  }, {"truth_pair_pt"})
        .Define("t_eta", [beta](float v) { return beta(v); }, {"truth_pair_eta"})
        .Define("t_m",   [bm](float v)   { return bm(v);   }, {"truth_minv"})
        .Define("r_pt",  [bpt](float v)  { return bpt(v);  }, {"pair_pt"})
        .Define("r_eta", [beta](float v) { return beta(v); }, {"pair_eta"})
        .Define("r_m",   [bm](float v)   { return bm(v);   }, {"minv"})
        .Define("inT", "pass_signal_truth && t_pt >= 0 && t_eta >= 0 && t_m >= 0")
        .Define("t_cell", Form("(t_pt >= 0 && t_eta >= 0 && t_m >= 0) ? (t_pt * %d + t_eta) * %d + t_m + 0.5 : -0.5", neta, nm))
        .Define("r_cell", Form("(r_pt >= 0 && r_eta >= 0 && r_m >= 0) ? (r_pt * %d + r_eta) * %d + r_m + 0.5 : -0.5", neta, nm))
        .Define("w_sf_tight",  "weight * m1.eff_sf_tight * m2.eff_sf_tight")
        .Define("w_sf_medium", "weight * m1.eff_sf_medium * m2.eff_sf_medium");

    const ROOT::RDF::TH1DModel mcell("", "", ncell, 0.0, ncell);
    const ROOT::RDF::TH2DModel mcell2("", "", ncell, 0.0, ncell, ncell, 0.0, ncell);
    auto M1 = [&](const std::string& n) { auto m = mcell;  m.fName = n; return m; };
    auto M2 = [&](const std::string& n) { auto m = mcell2; m.fName = n; return m; };
    const ROOT::RDF::TH1DModel mpt("", ";p_{T}^{pair} [GeV]", npt, ept.data());
    auto MP = [&](const std::string& n) { auto m = mpt; m.fName = n; return m; };

    std::vector<ROOT::RDF::RResultPtr<TH1D>> h1;
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h2;
    auto dT = df.Filter("inT");
    h1.push_back(dT.Histo1D(M1("h_cell_minv_T"), "t_cell", "weight"));
    h1.push_back(dT.Histo1D(MP("h_pair_pt_T"), "truth_pair_pt", "weight"));
    std::vector<ROOT::RDF::RResultPtr<double>> sf_sum, w_sum;
    for (const std::string wp : {"tight", "medium"}) {
        auto dR  = df.Filter("pair_pass_" + wp + " && pass_signal_reco");
        auto dTR = dR.Filter("inT");
        const std::string wsf = "w_sf_" + wp;
        h1.push_back(dR.Histo1D(M1("h_cell_minv_R_" + wp), "r_cell", "weight"));
        h1.push_back(dTR.Histo1D(M1("h_cell_minv_TR_" + wp), "t_cell", "weight"));
        h2.push_back(dTR.Histo2D(M2("h_cell_minv_M_" + wp), "t_cell", "r_cell", "weight"));
        h1.push_back(dR.Histo1D(M1("h_cell_minv_R_sf_" + wp), "r_cell", wsf));
        h1.push_back(dTR.Histo1D(M1("h_cell_minv_TR_sf_" + wp), "t_cell", wsf));
        h2.push_back(dTR.Histo2D(M2("h_cell_minv_M_sf_" + wp), "t_cell", "r_cell", wsf));
        h1.push_back(df.Filter("pair_pass_" + wp + " && pass_signal_reco_nogap")
                       .Histo1D(MP("h_pair_pt_Rnogap_" + wp), "pair_pt", "weight"));
        sf_sum.push_back(dR.Sum<double>(wsf));
        w_sum.push_back(dR.Sum<double>("weight"));
    }
    // SF fill sanity: muons of an R pair pass the WP, so their SF must be filled (> 0, never the default 1
    // for every muon). Count R muons whose SF is exactly 1 (unfilled default).
    auto n_sf1_tight = df.Filter("pair_pass_tight && pass_signal_reco")
                         .Filter("m1.eff_sf_tight == 1.f || m2.eff_sf_tight == 1.f").Count();
    auto n_R_tight = df.Filter("pair_pass_tight && pass_signal_reco").Count();

    TFile fout(out_file.c_str(), "RECREATE");
    for (auto& h : h1) h->Write();
    for (auto& h : h2) h->Write();
    TNamed("definition",
           "Cells as pair_bin_migration_pp24_full.root (cell_minv, (pT,eta^pair,m)); _sf_<wp> = weight*SF_1*SF_2 "
           "(m_i.eff_sf_<wp>) on the R leg; h_pair_pt_Rnogap_<wp> = R without the single-muon q*eta gap windows "
           "(|eta^pair|<2.2 kept); docs/tracking/crossx_raa_reco_binned_propagation_20261013.md").Write();
    TNamed("input", in_file.c_str()).Write();
    fout.Close();

    const char* wps[2] = {"tight", "medium"};
    for (int i = 0; i < 2; ++i)
        std::cout << Form("[sf] %s: <SF_1 SF_2> over R = %.5f\n", wps[i], *sf_sum[i] / *w_sum[i]);
    std::cout << "[sf] Tight R pairs with an unfilled (==1) leg SF: " << *n_sf1_tight << " / " << *n_R_tight << "\n";

    // ---- set identity with the nominal-tree migration file ----
    std::unique_ptr<TFile> fm(TFile::Open(mig_file.c_str(), "READ"));
    std::unique_ptr<TFile> fo(TFile::Open(out_file.c_str(), "READ"));
    if (!fm || fm->IsZombie()) throw std::runtime_error("cannot open " + mig_file);
    auto mine = [&](const std::string& n) {
        auto* h = dynamic_cast<TH1*>(fo->Get(n.c_str()));
        if (!h) throw std::runtime_error("missing " + n);
        return h;
    };
    CompareIdentical(mine("h_cell_minv_T"), fm.get(), "h_cell_minv_T");
    for (const std::string wp : {"tight", "medium"}) {
        CompareIdentical(mine("h_cell_minv_R_" + wp), fm.get(), "h_cell_minv_R_" + wp);
        CompareIdentical(mine("h_cell_minv_TR_" + wp), fm.get(), "h_cell_minv_TR_" + wp);
        CompareIdentical(mine("h_cell_minv_M_" + wp), fm.get(), "h_cell_minv_M_" + wp);
    }
    std::cout << "[done] " << out_file << "\n";
}
