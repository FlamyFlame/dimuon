// =================================================================================================
// fill_pbpb_muon_pt45_diag_counts.cxx
//
// TEMPORARY DIAGNOSTIC — muon reconstructed pT > 4.5 GeV cut study.
// docs/tracking/muon_pt45_cut_diagnostic.md. Human decision pending on raising the muon
// reconstructed-pT cut from the current 4 GeV (NTupleProcessingCode/DimuonDataAlgCoreT.c:596) to
// 4.5 GeV. If adopted, the permanent change belongs in NTuple processing (this file's cuts must
// then be deleted along with the matching block in RDFBasedHistFillingPP.cxx), every data and MC
// result must be rerun, and THIS FILE must be DELETED.
//
// WHY STANDALONE, not a block inside RDFBasedHistFillingPbPb::FillHistogramsCrossx (as done for
// pp24 in RDFBasedHistFillingPP.cxx): as of 2026-09, PbPb's trigger-weighted crossx event loop
// THROWS for some muons -- EvaluateSingleMuonEffcyPtFitted has no fitted turn-on, because PbPb's
// signal region (m1.charge*m1.eta < 2.2 below) has NOT migrated to the pp24 fiducial-gap cut that
// keeps every muon inside the fitted coarse q*eta range (pp24_crossx_rerun_2026_08.md). This is a
// PRE-EXISTING, unrelated, already-tracked bug: docs/tracking/pp_trig_eff_highpt_jump.md, ACTIVE,
// BLOCKED ON A USER DECISION. ROOT's RDF RunGraphs shares ONE event loop across every result
// booked from the same source tree, so that throw -- silently swallowed by ROOT, per
// reference_root_swallows_rdf_exceptions -- kills EVERY histogram tied to that loop, not just the
// trigger-weighted ones, which is why an in-class additive block (tried first, reverted) could not
// work no matter which pre-existing node it derived from.
//
// This macro therefore reads the ntuple-processing OUTPUT (muon_pair_tree_sign2 = OS, per
// reference_sign_convention) directly via a fresh, independent TChain/RDataFrame -- entirely
// decoupled from RDFBasedHistFillingPbPb and its poisoned trigger-efficiency machinery -- and
// mirrors the EXACT current PbPb signal selection (RDFBasedHistFillingPbPb.cxx:977,987) plus the
// Tight-WP filter (:987, "pair_pass_tight"). This is a pure raw-COUNT study: no trigger/reco
// efficiency, no luminosity scaling, no AMI/isospin weight -- none of that machinery is touched or
// needed. Combining the 3 years is a plain TChain union of the OS trees: for an UNWEIGHTED count,
// that is exactly the "simple sum, not luminosity-weighted" rule
// SingleBCrossxPlotterPbPbCombined::GetHistObject already uses for "_counts" histograms.
//
// TODO (not this diagnostic's scope): the PbPb crossx trigger-weighted hist filling is currently
// BROKEN for the reason above -- see docs/tracking/pp_trig_eff_highpt_jump.md (blocked on a user
// decision) and docs/tracking/muon_pt45_cut_diagnostic.md Remaining Work.
//
// Output: /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2/
//           histograms_pbpb_23_24_25_muon_pt45_diag_counts.root
//   h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts_diag_baseline  (muon pT > 4 GeV, current cut)
//   h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts_diag_mupt45    (muon pT > 4.5 GeV, candidate)
//   axes identical to pp24's h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts (ParamsSet::pT_bins_120
//   x ParamsSet::N_PAIR_ETA_CROSSX_BINS) -- .claude/CLAUDE.md Binnings: no new binning invented.
//
// Usage: root -l -b -q 'fill_pbpb_muon_pt45_diag_counts.cxx+()'
// =================================================================================================

#include <stdexcept>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TH2D.h"
#include "TSystem.h"
#include "ROOT/RDataFrame.hxx"

#include "../MuonObjectsParamsAndHelpers/ParamsSet.h"

void fill_pbpb_muon_pt45_diag_counts()
{
    const std::vector<int> years = {23, 24, 25};
    const std::string data_base = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data";

    TChain chain("muon_pair_tree_sign2");  // OS (reference_sign_convention.md)
    for (int yr : years) {
        const std::string path = data_base + "/pbpb_20" + std::to_string(yr)
            + "/muon_pairs_pbpb_20" + std::to_string(yr) + "_single_mu4_mindR_0_02.root";
        if (gSystem->AccessPathName(path.c_str())) {
            throw std::runtime_error("fill_pbpb_muon_pt45_diag_counts: input not found: " + path);
        }
        chain.Add(path.c_str());
        std::cout << "[INFO] Added " << path << std::endl;
    }

    ROOT::RDataFrame df(chain);

    // EXACT current PbPb signal region (RDFBasedHistFillingPbPb.cxx:977) + Tight WP (:987) --
    // mirrored verbatim, not re-derived. Deliberately NOT the pp24 fiducial-gap cut: pp and PbPb
    // currently use DIFFERENT signal regions (signal_selection_change_impact.md Sec.0) and this
    // diagnostic does not harmonize them.
    const std::string signal_cuts =
        "pair_pass_tight && minv > 1.08 && minv < 2.9 && pair_pt > 8"
        " && m1.charge*m1.eta < 2.2 && m2.charge*m2.eta < 2.2";

    ROOT::RDF::RNode df_signal = df.Filter(signal_cuts, "pbpb_signal_region");
    ROOT::RDF::RNode df_baseline = df_signal.Filter("m1.pt > 4 && m2.pt > 4", "diag_baseline_mupt4");
    ROOT::RDF::RNode df_mupt45   = df_signal.Filter("m1.pt > 4.5 && m2.pt > 4.5", "diag_mupt45");

    ParamsSet pms;
    const int npt = (int)(pms.pT_bins_120.size() - 1);
    const double* ptbins = pms.pT_bins_120.data();

    auto h_baseline = df_baseline.Histo2D(
        ROOT::RDF::TH2DModel("h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts_diag_baseline",
            ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins,
            ParamsSet::N_PAIR_ETA_CROSSX_BINS, ParamsSet::PAIR_ETA_CROSSX_MIN, ParamsSet::PAIR_ETA_CROSSX_MAX),
        "pair_pt", "pair_eta");
    auto h_mupt45 = df_mupt45.Histo2D(
        ROOT::RDF::TH2DModel("h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts_diag_mupt45",
            ";p_{T}^{pair} [GeV];#eta^{pair}", npt, ptbins,
            ParamsSet::N_PAIR_ETA_CROSSX_BINS, ParamsSet::PAIR_ETA_CROSSX_MIN, ParamsSet::PAIR_ETA_CROSSX_MAX),
        "pair_pt", "pair_eta");

    const double n_baseline = h_baseline->Integral(0, -1, 0, -1);
    const double n_mupt45   = h_mupt45->Integral(0, -1, 0, -1);
    std::cout << "[INFO] PbPb 23+24+25 combined, signal region, Tight WP:" << std::endl;
    std::cout << "[INFO]   muon pT>4.0 GeV (current) count = " << n_baseline << std::endl;
    std::cout << "[INFO]   muon pT>4.5 GeV (candidate) count = " << n_mupt45 << std::endl;
    if (n_baseline > 0) {
        std::cout << "[INFO]   pct decrease = " << 100.0 * (n_baseline - n_mupt45) / n_baseline
                  << " %" << std::endl;
    }

    const std::string out_dir = data_base + "/pbpb_run2";
    gSystem->mkdir(out_dir.c_str(), true);
    const std::string out_path = out_dir + "/histograms_pbpb_23_24_25_muon_pt45_diag_counts.root";
    TFile fout(out_path.c_str(), "RECREATE");
    h_baseline->Write();
    h_mupt45->Write();
    fout.Close();
    std::cout << "[INFO] Saved: " << out_path << std::endl;
}
