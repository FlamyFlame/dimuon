// =================================================================================================
// fill_pp24_fullsim_pair_bin_migration.C
//
// Bin-migration ingredients for the bin-by-bin vs full-unfolding decision
// (docs/tracking/pair_bin_migration_bin_by_bin_feasibility.md, Physics Procedure §2-§3a).
//
// Input : the ntuple-processing output muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_full.root,
//         tree muon_pair_tree_sign2 (truth OS pairs), weight `weight` -- exactly the RDF class's
//         `df_op_weighted` (RDFBasedHistFillingBaseClass::CreateBaseRDFsBaseCommon, no extra filter).
// Cuts  : the SAME expression builders as RDFBasedHistFillingPythiaFullsim.cxx
//         (ParamsSet::SignalMinvCutExpr / SignalPairPtCutExpr, RDFBasedHistFillingPythiaFullsim::RecoGapCutExpr).
// Axes  : ParamsSet::pair_pt_coarse_bins, CommonEffcyConfig::pair_eta_proj_ranges_coarse_incl_gap,
//         ParamsSet::SignalMinvRecoEffcyBins(), CommonEffcyConfig::dr_bins_edges_for_reco_effcy.
//
//   T = truth signal region (mass window, pair pT) with truth (pT, eta^pair, m) inside the axes
//   R = pair passes WP + RECO signal selection incl. gap cuts (= the applied eps_reco numerator's
//       reco leg); a reco value outside an axis goes to the overflow = migrated out of range
//
// Written per variable x (1D: pair_pt, pair_eta, minv, dr; 3D cell index: cell_minv = (pT,eta,m),
// cell_dr = (pT,eta,dR)):
//   h_<x>_T            N[T]        vs truth x           (WP-independent)
//   h_<x>_TR_<wp>      N[T and R]  vs truth x
//   h_<x>_M_<wp>       N[T and R]  truth x (X) vs reco x (Y)
//   h_<x>_R_<wp>       N[R]        vs reco x  (truth anywhere)
//   h_<x>_Rout_<wp>    N[R and !T] vs reco x  (feed-in from outside the measured range)
// 3D cells are flattened as (i_pt * N_eta + i_eta) * N_z + i_z, 0-based.
//
// Closure: N[T] and N[T and R] in cell_dr / cell_minv are compared with the RDF-class 3D
// `h_truth_{dr,minv}_effcy_vs_truth_pair_eta_coarse_vs_truth_pair_pt_coarse_op{_pass_signal_truth,
// _pass_<wp>_and_signal_truth_and_reco}` cell by cell; the macro THROWS on a mismatch.
//
// Usage (from Analysis/plotting_codes/reco_effcy/bin_migration/):
//   root -l -b -q 'fill_pp24_fullsim_pair_bin_migration.C+'
// =================================================================================================

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>
#include <TAxis.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TSystem.h>

#include "../../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../../RDFBasedHistFilling/RDFBasedHistFillingPythia.h"

namespace {

struct Axes {
    std::vector<double> pt, eta, minv, dr;
};

Axes CanonicalAxes()
{
    static const ParamsSet pms{};
    static const CommonEffcyConfig cfg{};
    Axes a;
    a.pt = pms.pair_pt_coarse_bins;
    a.eta.push_back(cfg.pair_eta_proj_ranges_coarse_incl_gap.front().first);
    for (const auto& r : cfg.pair_eta_proj_ranges_coarse_incl_gap) a.eta.push_back(r.second);
    a.minv = ParamsSet::SignalMinvRecoEffcyBins();
    a.dr.assign(cfg.dr_bins_edges_for_reco_effcy.begin(), cfg.dr_bins_edges_for_reco_effcy.end());
    return a;
}

// 0-based bin index on the given edges (TAxis convention: low edge inclusive), -1 outside.
struct BinIdx {
    TAxis ax;
    explicit BinIdx(const std::vector<double>& e) : ax((int)e.size() - 1, e.data()) {}
    int operator()(double v) const
    {
        const int b = ax.FindFixBin(v);
        return (b >= 1 && b <= ax.GetNbins()) ? b - 1 : -1;
    }
};

void CompareCells(TFile* fh, const TH1D* mine, const std::string& name3d, const std::string& what)
{
    auto* h = dynamic_cast<TH3D*>(fh->Get(name3d.c_str()));
    if (!h) throw std::runtime_error("closure: missing " + name3d);
    const int nx = h->GetNbinsX(), ny = h->GetNbinsY(), nz = h->GetNbinsZ();
    if (nx * ny * nz != mine->GetNbinsX())
        throw std::runtime_error("closure: cell count differs for " + what);
    double max_rel = 0.0;
    for (int ix = 1; ix <= nx; ++ix)
        for (int iy = 1; iy <= ny; ++iy)
            for (int iz = 1; iz <= nz; ++iz) {
                const int k = ((ix - 1) * ny + (iy - 1)) * nz + (iz - 1);
                const double a = h->GetBinContent(ix, iy, iz), b = mine->GetBinContent(k + 1);
                if (a == 0.0 && b == 0.0) continue;
                max_rel = std::max(max_rel, std::fabs(a - b) / std::max(std::fabs(a), std::fabs(b)));
            }
    std::cout << "[closure] " << what << " vs " << name3d << ": max |rel diff| = " << max_rel << "\n";
    if (max_rel > 1e-6) throw std::runtime_error("closure FAILED for " + what);
}

} // namespace

void fill_pp24_fullsim_pair_bin_migration()
{
    const Axes A = CanonicalAxes();
    const int npt = A.pt.size() - 1, neta = A.eta.size() - 1, nm = A.minv.size() - 1,
              ndr = A.dr.size() - 1;
    const int ncell_m = npt * neta * nm, ncell_dr = npt * neta * ndr;

    const std::string sample_dir = FullSimSampleInputDir(FullSimSampleType::pp, /*is_test_sample=*/false);
    const std::string in_file  = sample_dir + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_full.root";
    const std::string hist_file = sample_dir + "histograms_pythia_fullsim_pp24_no_data_resonance_cuts_full.root";
    const std::string out_dir  = FullSimRecoEffDir(sample_dir);
    const std::string out_file = out_dir + "pair_bin_migration_pp24_full.root";
    gSystem->mkdir(out_dir.c_str(), true);

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df0("muon_pair_tree_sign2", in_file);

    const BinIdx bpt(A.pt), beta(A.eta), bm(A.minv), bdr(A.dr);
    // Signal definitions: the RDF class's own strings (RDFBasedHistFillingPythiaFullsim.cxx).
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
        .Define("t_dr",   [bdr](float v)  { return bdr(v);  }, {"truth_dr"})
        .Define("r_pt",   [bpt](float v)  { return bpt(v);  }, {"pair_pt"})
        .Define("r_eta",  [beta](float v) { return beta(v); }, {"pair_eta"})
        .Define("r_m",    [bm](float v)   { return bm(v);   }, {"minv"})
        .Define("r_dr",   [bdr](float v)  { return bdr(v);  }, {"dr"})
        .Define("inT", "pass_signal_truth && t_pt >= 0 && t_eta >= 0 && t_m >= 0")
        .Define("t_cell_minv", Form("(t_pt >= 0 && t_eta >= 0 && t_m >= 0) ? (t_pt * %d + t_eta) * %d + t_m : -1", neta, nm))
        .Define("r_cell_minv", Form("(r_pt >= 0 && r_eta >= 0 && r_m >= 0) ? (r_pt * %d + r_eta) * %d + r_m : -1", neta, nm))
        .Define("t_cell_dr",   Form("(t_pt >= 0 && t_eta >= 0 && t_dr >= 0) ? (t_pt * %d + t_eta) * %d + t_dr : -1", neta, ndr))
        .Define("r_cell_dr",   Form("(r_pt >= 0 && r_eta >= 0 && r_dr >= 0) ? (r_pt * %d + r_eta) * %d + r_dr : -1", neta, ndr))
        // index columns as doubles (+0.5 -> bin centre) for the cell histograms
        .Define("t_cell_minv_d", "t_cell_minv + 0.5").Define("r_cell_minv_d", "r_cell_minv + 0.5")
        .Define("t_cell_dr_d",   "t_cell_dr + 0.5").Define("r_cell_dr_d",   "r_cell_dr + 0.5");

    struct Var {
        std::string key, tcol, rcol;
        std::vector<double> edges;  // empty -> cell index axis with ncell bins
        int ncell;
    };
    const std::vector<Var> vars = {
        {"pair_pt",   "truth_pair_pt",  "pair_pt",  A.pt,   0},
        {"pair_eta",  "truth_pair_eta", "pair_eta", A.eta,  0},
        {"minv",      "truth_minv",     "minv",     A.minv, 0},
        {"dr",        "truth_dr",       "dr",       A.dr,   0},
        {"cell_minv", "t_cell_minv_d",  "r_cell_minv_d", {}, ncell_m},
        {"cell_dr",   "t_cell_dr_d",    "r_cell_dr_d",   {}, ncell_dr},
    };
    auto model1 = [](const Var& v, const std::string& n) {
        return v.edges.empty() ? ROOT::RDF::TH1DModel(n.c_str(), "", v.ncell, 0.0, v.ncell)
                               : ROOT::RDF::TH1DModel(n.c_str(), "", (int)v.edges.size() - 1, v.edges.data());
    };
    auto model2 = [](const Var& v, const std::string& n) {
        return v.edges.empty()
            ? ROOT::RDF::TH2DModel(n.c_str(), "", v.ncell, 0.0, v.ncell, v.ncell, 0.0, v.ncell)
            : ROOT::RDF::TH2DModel(n.c_str(), "", (int)v.edges.size() - 1, v.edges.data(),
                                   (int)v.edges.size() - 1, v.edges.data());
    };

    std::vector<ROOT::RDF::RResultPtr<TH1D>> h1;
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h2;
    auto dT = df.Filter("inT");
    for (const auto& v : vars)
        h1.push_back(dT.Histo1D(model1(v, "h_" + v.key + "_T"), v.tcol, "weight"));
    for (const std::string wp : {"tight", "medium"}) {
        // R = the APPLIED eps_reco reco selection, nothing more: a reco value outside an axis
        // (reco pair pT > top edge; reco dR > top edge is possible in principle) lands in the overflow (cell index -1 -> underflow)
        // and counts as migrated OUT of the measured range.
        auto dR  = df.Filter("pair_pass_" + wp + " && pass_signal_reco");
        auto dTR = dR.Filter("inT");
        auto dRo = dR.Filter("!inT");
        for (const auto& v : vars) {
            h1.push_back(dTR.Histo1D(model1(v, "h_" + v.key + "_TR_" + wp), v.tcol, "weight"));
            h2.push_back(dTR.Histo2D(model2(v, "h_" + v.key + "_M_" + wp), v.tcol, v.rcol, "weight"));
            h1.push_back(dR.Histo1D(model1(v, "h_" + v.key + "_R_" + wp), v.rcol, "weight"));
            h1.push_back(dRo.Histo1D(model1(v, "h_" + v.key + "_Rout_" + wp), v.rcol, "weight"));
        }
    }
    // Bookkeeping: T pairs outside the dR axis (must be 0: dR is kinematically bounded) and
    // Tight-selected reco pairs outside the (pT, eta, m) axes.
    auto n_T_dr_out  = dT.Filter("t_dr < 0").Sum<double>("weight");
    auto n_R_axisout = df.Filter("pair_pass_tight && pass_signal_reco && !(r_pt >= 0 && r_eta >= 0 && r_m >= 0)").Sum<double>("weight");
    auto n_R_dr_out  = df.Filter("pair_pass_tight && pass_signal_reco && r_dr < 0").Sum<double>("weight");
    auto n_R_sel     = df.Filter("pair_pass_tight && pass_signal_reco").Sum<double>("weight");

    TFile fout(out_file.c_str(), "RECREATE");
    for (auto& h : h1) h->Write();   // RDF weighted fills already carry Sumw2
    for (auto& h : h2) h->Write();
    // Axes and the definition, for the plotter (it re-verifies them against the canonical source).
    TNamed("definition",
           "T = truth signal region (truth mass window && truth pair pT) && truth (pT, eta^pair, m) inside "
           "the canonical axes; R = pair_pass_<wp> && pass_signal_reco (reco mass window, pair pT, gap "
           "cuts incl. |eta^pair|) = exactly the applied eps_reco reco leg; a reco value outside an axis "
           "(only reco pT > top edge occurs) is in the overflow (underflow for the cell-index histograms) = migrated out; category = all truth OS "
           "pairs (_op); docs/tracking/pair_bin_migration_bin_by_bin_feasibility.md").Write();
    TNamed("input", in_file.c_str()).Write();
    fout.Close();

    std::cout << "[info] T pairs outside the dR axis (weight): " << *n_T_dr_out << "\n"
              << "[info] Tight R pairs outside the dR axis: " << 100.0 * *n_R_dr_out / *n_R_sel << " %\n"
              << "[info] Tight R pairs outside the (pT, eta, m) axes: " << 100.0 * *n_R_axisout / *n_R_sel << " %\n";
    if (*n_T_dr_out != 0.0) throw std::runtime_error("T pairs outside the dR axis -- dR axis does not cover T");

    // ---- closure against the RDF-class 3D eps_reco legs ----
    TFile fh(hist_file.c_str(), "READ");
    TFile fm(out_file.c_str(), "READ");
    auto get = [&](const std::string& n) {
        auto* h = dynamic_cast<TH1D*>(fm.Get(n.c_str()));
        if (!h) throw std::runtime_error("missing " + n);
        return h;
    };
    const std::string b3 = "_vs_truth_pair_eta_coarse_vs_truth_pair_pt_coarse_op";
    CompareCells(&fh, get("h_cell_dr_T"),   "h_truth_dr_effcy" + b3 + "_pass_signal_truth",   "cell_dr T");
    CompareCells(&fh, get("h_cell_minv_T"), "h_truth_minv_effcy" + b3 + "_pass_signal_truth", "cell_minv T");
    for (const std::string wp : {"tight", "medium"}) {
        CompareCells(&fh, get("h_cell_dr_TR_" + wp),   "h_truth_dr_effcy" + b3 + "_pass_" + wp + "_and_signal_truth_and_reco",   "cell_dr TR " + wp);
        CompareCells(&fh, get("h_cell_minv_TR_" + wp), "h_truth_minv_effcy" + b3 + "_pass_" + wp + "_and_signal_truth_and_reco", "cell_minv TR " + wp);
    }
    std::cout << "[done] " << out_file << "\n";
}
