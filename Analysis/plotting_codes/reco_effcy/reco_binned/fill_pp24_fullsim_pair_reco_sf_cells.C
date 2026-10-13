// =================================================================================================
// fill_pp24_fullsim_pair_reco_sf_cells.C
//
// SF-weighted (pair pT, eta^pair, m_mumu) cell ingredients of the data-corrected reco-binned pair
// reco efficiency, and the NO-GAP reco leg used to scale the Pb+Pb Run-2 placeholder
// (docs/tracking/crossx_raa_reco_binned_propagation_20261013.md, Physics Procedure §2, §3a):
//
//   eps^{reco,SF}_j = sum_{R_j} w * SF_1 * SF_2 / N[T_j]       SF_i = m_i.eff_sf_<wp>
//
// SF below the tool's validity (D5, user 2026-10-13: "nearest valid SF"): the skim zeroes the SF when the
// MuonEfficiencyScaleFactors tool returns a non-Ok code (TrigRates.cxx:1322-1341) and the NTP then stores 1
// (PythiaFullSimExtras.c:358-361). This happens for EVERY muon with pT < 5 GeV and for no muon above.
// Such a muon gets the SF of the nearest valid pT slice, [5, 6) GeV, at its eta: the MC-weighted mean
// filled SF of WP-passing R-pair muons in that slice, in eta bins of kSfEtaBinWidth (a lookup granularity,
// not a result binning). SF = 1 is kept as the systematic variation (h_*_sf1_<wp> histograms).
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
#include <TProfile.h>
#include <map>
#include <TNamed.h>
#include <TROOT.h>
#include <TSystem.h>

#include "../../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../../RDFBasedHistFilling/RDFBasedHistFillingPythia.h"
#include "../../../Utilities/PairRecoEffRecoBinnedEvaluator.h"

namespace {

constexpr double kSfValidPtMin  = 5.0;   // GeV: lowest pT at which the skim's SF tool returns a valid SF
constexpr double kSfValidPtMax  = 6.0;   // GeV: upper edge of the "nearest valid" slice
constexpr double kSfEtaBinWidth = 0.1;   // eta granularity of the nearest-valid SF lookup
constexpr double kSfEtaMax      = 2.5;

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
        .Define("w_sf1_tight",  "weight * m1.eff_sf_tight * m2.eff_sf_tight")      // SF = 1 below 5 GeV (systematic)
        .Define("w_sf1_medium", "weight * m1.eff_sf_medium * m2.eff_sf_medium");

    // ---- nearest-valid SF lookup: <SF>(eta) of filled WP muons with 5 <= pT < 6 GeV (first event loop) ----
    const int neta_sf = (int)std::lround(2 * kSfEtaMax / kSfEtaBinWidth);
    std::map<std::string, std::shared_ptr<TProfile>> sf_lookup;
    for (const std::string wp : {"tight", "medium"}) {
        const std::string pass = "pair_pass_" + wp + " && pass_signal_reco";
        const std::string sel = Form("m%%d.pt >= %g && m%%d.pt < %g && m%%d.eff_sf_%s != 1.f", kSfValidPtMin, kSfValidPtMax, wp.c_str());
        auto p1 = df.Filter(pass).Filter(Form(sel.c_str(), 1, 1, 1))
                    .Profile1D({("p1_" + wp).c_str(), "", neta_sf, -kSfEtaMax, kSfEtaMax}, "m1.eta", "m1.eff_sf_" + wp, "weight");
        auto p2 = df.Filter(pass).Filter(Form(sel.c_str(), 2, 2, 2))
                    .Profile1D({("p2_" + wp).c_str(), "", neta_sf, -kSfEtaMax, kSfEtaMax}, "m2.eta", "m2.eff_sf_" + wp, "weight");
        auto prof = std::shared_ptr<TProfile>(static_cast<TProfile*>(p1->Clone(("sf_lookup_" + wp).c_str())));
        prof->SetDirectory(nullptr);
        prof->Add(p2.GetPtr());
        for (int b = 1; b <= neta_sf; ++b)
            if (!(prof->GetBinEntries(b) > 0) && std::fabs(prof->GetBinCenter(b)) < 2.4)
                throw std::runtime_error("nearest-valid SF lookup: empty eta bin " + std::to_string(b) + " (" + wp + ")");
        sf_lookup[wp] = prof;
    }
    auto sf_fix = [](const std::shared_ptr<TProfile>& prof) {
        return [prof](float pt, float eta, float sf) {
            if (sf != 1.f || pt >= kSfValidPtMin) return (double)sf;
            const int b = std::min(std::max(prof->FindFixBin(eta), 1), prof->GetNbinsX());
            if (!(prof->GetBinEntries(b) > 0))
                throw std::runtime_error("nearest-valid SF lookup: no filled [5,6) GeV muon at eta " + std::to_string(eta));
            return prof->GetBinContent(b);
        };
    };
    auto df_sf = df
        .Define("sf1_tight",  sf_fix(sf_lookup["tight"]),  {"m1.pt", "m1.eta", "m1.eff_sf_tight"})
        .Define("sf2_tight",  sf_fix(sf_lookup["tight"]),  {"m2.pt", "m2.eta", "m2.eff_sf_tight"})
        .Define("sf1_medium", sf_fix(sf_lookup["medium"]), {"m1.pt", "m1.eta", "m1.eff_sf_medium"})
        .Define("sf2_medium", sf_fix(sf_lookup["medium"]), {"m2.pt", "m2.eta", "m2.eff_sf_medium"})
        .Define("w_sf_tight",  "weight * sf1_tight * sf2_tight")
        .Define("w_sf_medium", "weight * sf1_medium * sf2_medium");

    const ROOT::RDF::TH1DModel mcell("", "", ncell, 0.0, ncell);
    const ROOT::RDF::TH2DModel mcell2("", "", ncell, 0.0, ncell, ncell, 0.0, ncell);
    auto M1 = [&](const std::string& n) { auto m = mcell;  m.fName = n; return m; };
    auto M2 = [&](const std::string& n) { auto m = mcell2; m.fName = n; return m; };
    const ROOT::RDF::TH1DModel mpt("", ";p_{T}^{pair} [GeV]", npt, ept.data());
    auto MP = [&](const std::string& n) { auto m = mpt; m.fName = n; return m; };

    std::vector<ROOT::RDF::RResultPtr<TH1D>> h1;
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h2;
    auto dT = df_sf.Filter("inT");
    h1.push_back(dT.Histo1D(M1("h_cell_minv_T"), "t_cell", "weight"));
    h1.push_back(dT.Histo1D(MP("h_pair_pt_T"), "truth_pair_pt", "weight"));
    std::vector<ROOT::RDF::RResultPtr<double>> sf_sum, w_sum;
    for (const std::string wp : {"tight", "medium"}) {
        auto dR  = df_sf.Filter("pair_pass_" + wp + " && pass_signal_reco");
        auto dTR = dR.Filter("inT");
        const std::string wsf = "w_sf_" + wp;
        h1.push_back(dR.Histo1D(M1("h_cell_minv_R_sf1_" + wp), "r_cell", "w_sf1_" + wp));   // SF = 1 below 5 GeV
        h1.push_back(dR.Histo1D(M1("h_cell_minv_R_" + wp), "r_cell", "weight"));
        h1.push_back(dTR.Histo1D(M1("h_cell_minv_TR_" + wp), "t_cell", "weight"));
        h2.push_back(dTR.Histo2D(M2("h_cell_minv_M_" + wp), "t_cell", "r_cell", "weight"));
        h1.push_back(dR.Histo1D(M1("h_cell_minv_R_sf_" + wp), "r_cell", wsf));
        h1.push_back(dTR.Histo1D(M1("h_cell_minv_TR_sf_" + wp), "t_cell", wsf));
        h2.push_back(dTR.Histo2D(M2("h_cell_minv_M_sf_" + wp), "t_cell", "r_cell", wsf));
        h1.push_back(df_sf.Filter("pair_pass_" + wp + " && pass_signal_reco_nogap")
                       .Histo1D(MP("h_pair_pt_Rnogap_" + wp), "pair_pt", "weight"));
        sf_sum.push_back(dR.Sum<double>(wsf));
        w_sum.push_back(dR.Sum<double>("weight"));
    }
    // SF fill sanity: an unfilled SF (default 1) must occur ONLY below 5 GeV -- else the nearest-valid
    // substitution would be hiding a different failure.
    auto n_sf1_tight = df.Filter("pair_pass_tight && pass_signal_reco")
                         .Filter("m1.eff_sf_tight == 1.f || m2.eff_sf_tight == 1.f").Count();
    auto n_sf1_above = df.Filter("pair_pass_tight && pass_signal_reco")
                         .Filter(Form("(m1.eff_sf_tight == 1.f && m1.pt >= %g) || (m2.eff_sf_tight == 1.f && m2.pt >= %g)",
                                      kSfValidPtMin, kSfValidPtMin)).Count();
    auto n_sf1_above_medium = df.Filter("pair_pass_medium && pass_signal_reco")
                         .Filter(Form("(m1.eff_sf_medium == 1.f && m1.pt >= %g) || (m2.eff_sf_medium == 1.f && m2.pt >= %g)",
                                      kSfValidPtMin, kSfValidPtMin)).Count();
    auto n_R_tight = df.Filter("pair_pass_tight && pass_signal_reco").Count();

    TFile fout(out_file.c_str(), "RECREATE");
    for (auto& h : h1) h->Write();
    for (auto& h : h2) h->Write();
    for (auto& kv : sf_lookup) kv.second->Write();
    TNamed("definition",
           "Cells as pair_bin_migration_pp24_full.root (cell_minv, (pT,eta^pair,m)); _sf_<wp> = weight*SF_1*SF_2 "
           "(m_i.eff_sf_<wp>; pT < 5 GeV -> nearest valid [5,6) GeV SF at the muon eta, sf_lookup_<wp>; _sf1_ = SF 1 "
           "there) on the R leg; h_pair_pt_Rnogap_<wp> = R without the single-muon q*eta gap windows "
           "(|eta^pair|<2.2 kept); docs/tracking/crossx_raa_reco_binned_propagation_20261013.md").Write();
    TNamed("input", in_file.c_str()).Write();
    fout.Close();

    const char* wps[2] = {"tight", "medium"};
    for (int i = 0; i < 2; ++i)
        std::cout << Form("[sf] %s: <SF_1 SF_2> over R = %.5f\n", wps[i], *sf_sum[i] / *w_sum[i]);
    std::cout << "[sf] Tight R pairs with an unfilled (==1) leg SF: " << *n_sf1_tight << " / " << *n_R_tight
              << "; of which with the unfilled leg at pT >= " << kSfValidPtMin << " GeV: " << *n_sf1_above << "\n";
    std::cout << "[sf] Medium R pairs with an unfilled leg at pT >= " << kSfValidPtMin << " GeV: " << *n_sf1_above_medium << "\n";
    if (*n_sf1_above != 0 || *n_sf1_above_medium != 0) throw std::runtime_error("unfilled SF above the tool validity edge -- investigate before using");

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
