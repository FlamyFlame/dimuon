// =================================================================================================
// fill_pp24_fullsim_pair_eta_smoothness_closure.C
//
// MC closure of the eta^pair smoothness test (docs/tracking/pair_reco_eff_bin_by_bin_reco_binned.md
// §3c): MC R (reco eta^pair) corrected by 1/eps^reco(reco pT, eta^pair, m) on the COARSE axes vs MC T
// (truth eta^pair), both in the FINE plot-only binning. Summed over each coarse cell they agree by
// construction (closure check, printed); inside a coarse cell the fine-bin structure shows the
// coarse-eta averaging effect, free of data statistics.
//
// Input / selection: EXACTLY the migration filler (fill_pp24_fullsim_pair_bin_migration.C) -- the
// ntuple-processing output muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_full.root, tree
// muon_pair_tree_sign2 (all truth OS pairs), weight `weight`, pass_signal_truth / pass_signal_reco
// from the SAME ParamsSet / RDFBasedHistFillingPythiaFullsim builders, T = inT, R = pair_pass_<wp> &&
// pass_signal_reco. The definitions below are copied verbatim from that filler.
// Also re-fills N[T_j], N[R_j] per (pT, eta^pair, m) cell and compares them with the migration file
// cell by cell (consistency of the migration file with the CURRENT ntuple); THROWS on a mismatch.
//
// Output: pair_eta_smoothness_common.h SmoothnessMCFile():
//   h_mcT_fine[_ptK], h_mcR_fine_raw_<wp>[_ptK], h_mcR_fine_corr_<wp>[_ptK]  (K = coarse pT bin;
//   T panels by TRUTH pT, R panels by RECO pT).
//
// Usage (from Analysis/plotting_codes/reco_effcy/reco_binned/):
//   root -l -b -q 'fill_pp24_fullsim_pair_eta_smoothness_closure.C+'
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
#include <TNamed.h>
#include <TROOT.h>
#include <TSystem.h>

#include "../../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../../RDFBasedHistFilling/RDFBasedHistFillingPythia.h"
#include "pair_eta_smoothness_common.h"

namespace {
struct BinIdx {   // == the migration filler's BinIdx
    TAxis ax;
    explicit BinIdx(const std::vector<double>& e) : ax((int)e.size() - 1, e.data()) {}
    int operator()(double v) const
    {
        const int b = ax.FindFixBin(v);
        return (b >= 1 && b <= ax.GetNbins()) ? b - 1 : -1;
    }
};
PairRecoEffRecoBinnedEvaluator* s_ev_tight = nullptr;
PairRecoEffRecoBinnedEvaluator* s_ev_medium = nullptr;
}  // namespace

void fill_pp24_fullsim_pair_eta_smoothness_closure()
{
    using EV = PairRecoEffRecoBinnedEvaluator;
    const std::vector<double> ept = EV::PtEdges(), eeta = EV::EtaEdges(), em = EV::MinvEdges();
    const int npt = ept.size() - 1, neta = eeta.size() - 1, nm = em.size() - 1, ncell = npt * neta * nm;
    const std::vector<double> fine = PairEtaFineSmoothnessBins();

    const std::string sample_dir = FullSimSampleInputDir(FullSimSampleType::pp, /*is_test_sample=*/false);
    const std::string in_file  = sample_dir + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_full.root";
    const std::string mig_file = FullSimRecoEffDir(sample_dir) + "pair_bin_migration_pp24_full.root";
    const std::string out_file = SmoothnessMCFile();
    gSystem->mkdir(SmoothnessDataDir().c_str(), true);

    s_ev_tight = new EV();  s_ev_tight->Load(RecoBinnedProductFile(), true);
    s_ev_medium = new EV(); s_ev_medium->Load(RecoBinnedProductFile(), false);

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df0("muon_pair_tree_sign2", in_file);
    const BinIdx bpt(ept), beta(eeta), bm(em);
    // ---- verbatim from fill_pp24_fullsim_pair_bin_migration.C ----
    const std::string sig_truth = ParamsSet::SignalMinvCutExpr("truth_minv") + " && "
                                + ParamsSet::SignalPairPtCutExpr("truth_pair_pt");
    const std::string sig_reco  = "(m1.reco_match && m2.reco_match) ? (" + ParamsSet::SignalMinvCutExpr("minv")
                                + " && " + ParamsSet::SignalPairPtCutExpr("pair_pt") + " && "
                                + RDFBasedHistFillingPythiaFullsim::RecoGapCutExpr() + ") : false";
    auto df = df0
        .Define("pass_signal_truth", sig_truth)
        .Define("pass_signal_reco",  sig_reco)
        .Define("t_pt",   [bpt](float v)  { return bpt(v);  }, {"truth_pair_pt"})
        .Define("t_eta",  [beta](float v) { return beta(v); }, {"truth_pair_eta"})
        .Define("t_m",    [bm](float v)   { return bm(v);   }, {"truth_minv"})
        .Define("r_pt",   [bpt](float v)  { return bpt(v);  }, {"pair_pt"})
        .Define("r_eta",  [beta](float v) { return beta(v); }, {"pair_eta"})
        .Define("r_m",    [bm](float v)   { return bm(v);   }, {"minv"})
        .Define("inT", "pass_signal_truth && t_pt >= 0 && t_eta >= 0 && t_m >= 0")
        .Define("t_cell_minv", Form("(t_pt >= 0 && t_eta >= 0 && t_m >= 0) ? (t_pt * %d + t_eta) * %d + t_m : -1", neta, nm))
        .Define("r_cell_minv", Form("(r_pt >= 0 && r_eta >= 0 && r_m >= 0) ? (r_pt * %d + r_eta) * %d + r_m : -1", neta, nm))
        .Define("t_cell_minv_d", "t_cell_minv + 0.5").Define("r_cell_minv_d", "r_cell_minv + 0.5")
    // ---- end verbatim ----
        .Define("w_corr_tight", [](float pt, float eta, float m, double w) {
            const double e = s_ev_tight->Eval(pt, eta, m); return e > 0 ? w / e : w; }, {"pair_pt", "pair_eta", "minv", "weight"})
        .Define("w_corr_medium", [](float pt, float eta, float m, double w) {
            const double e = s_ev_medium->Eval(pt, eta, m); return e > 0 ? w / e : w; }, {"pair_pt", "pair_eta", "minv", "weight"});

    std::vector<ROOT::RDF::RResultPtr<TH1D>> hs;
    auto fineModel = [&](const std::string& n) {
        return ROOT::RDF::TH1DModel(n.c_str(), ";#eta^{pair};", (int)fine.size() - 1, fine.data());
    };
    auto cellModel = [&](const std::string& n) { return ROOT::RDF::TH1DModel(n.c_str(), "", ncell, 0.0, ncell); };

    auto dT = df.Filter("inT");
    hs.push_back(dT.Histo1D(fineModel("h_mcT_fine"), "truth_pair_eta", "weight"));
    hs.push_back(dT.Histo1D(cellModel("h_chk_cell_T"), "t_cell_minv_d", "weight"));
    for (int k = 0; k < npt; ++k)
        hs.push_back(dT.Filter(Form("t_pt == %d", k)).Histo1D(fineModel(Form("h_mcT_fine_pt%d", k)), "truth_pair_eta", "weight"));
    for (const std::string wp : {"tight", "medium"}) {
        auto dR = df.Filter("pair_pass_" + wp + " && pass_signal_reco");
        const std::string wc = "w_corr_" + wp;
        hs.push_back(dR.Histo1D(fineModel("h_mcR_fine_raw_" + wp), "pair_eta", "weight"));
        hs.push_back(dR.Histo1D(fineModel("h_mcR_fine_corr_" + wp), "pair_eta", wc));
        hs.push_back(dR.Histo1D(cellModel("h_chk_cell_R_" + wp), "r_cell_minv_d", "weight"));
        for (int k = 0; k < npt; ++k) {
            auto d = dR.Filter(Form("r_pt == %d", k));
            hs.push_back(d.Histo1D(fineModel(Form("h_mcR_fine_raw_%s_pt%d", wp.c_str(), k)), "pair_eta", "weight"));
            hs.push_back(d.Histo1D(fineModel(Form("h_mcR_fine_corr_%s_pt%d", wp.c_str(), k)), "pair_eta", wc));
        }
    }

    TFile fout(out_file.c_str(), "RECREATE");
    for (auto& h : hs) h->Write();
    TNamed("input", in_file.c_str()).Write();
    TNamed("reco_binned_product", RecoBinnedProductFile().c_str()).Write();
    TNamed("definition", "MC closure: T = inT (truth eta^pair), R = pair_pass_<wp> && pass_signal_reco (reco eta^pair), "
                         "corr = weight / eps^reco(reco pT, eta^pair, m) coarse axes; selection verbatim from "
                         "fill_pp24_fullsim_pair_bin_migration.C; fine axis = ParamsSet::pair_eta_measured_bins").Write();
    fout.Close();
    s_ev_tight->PrintStats();
    s_ev_medium->PrintStats();

    // ---- consistency of the migration file with the current ntuple ----
    std::unique_ptr<TFile> fm(TFile::Open(mig_file.c_str(), "READ"));
    std::unique_ptr<TFile> fo(TFile::Open(out_file.c_str(), "READ"));
    auto cmp = [&](const std::string& mine, const std::string& theirs) {
        auto* a = (TH1D*)fo->Get(mine.c_str());
        auto* b = (TH1D*)fm->Get(theirs.c_str());
        if (!a || !b) throw std::runtime_error("consistency: missing " + mine + " / " + theirs);
        double mx = 0;
        for (int i = 0; i <= ncell + 1; ++i) {
            const double x = a->GetBinContent(i), y = b->GetBinContent(i);
            if (x == 0 && y == 0) continue;
            mx = std::max(mx, std::fabs(x - y) / std::max(std::fabs(x), std::fabs(y)));
        }
        std::cout << "[consistency] " << mine << " vs migration " << theirs << ": max |rel diff| = " << mx << "\n";
        if (mx > 1e-6) throw std::runtime_error("migration file NOT consistent with the current ntuple: " + theirs);
    };
    cmp("h_chk_cell_T", "h_cell_minv_T");
    cmp("h_chk_cell_R_tight", "h_cell_minv_R_tight");
    cmp("h_chk_cell_R_medium", "h_cell_minv_R_medium");
    std::cout << "[done] " << out_file << "\n";
}
