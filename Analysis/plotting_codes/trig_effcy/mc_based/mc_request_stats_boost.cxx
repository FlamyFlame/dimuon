// Projected statistics boost of the additional pp24 Pythia MC request in the 3 highest pair-pT x
// 3 |eta^pair| cells (docs/tracking/additional_mc_request_jo_test.md, Physics Procedure §3d, D5).
//
// REFERENCE: the existing pp24 fullsim FULL sample, Step-3 MC trigger-efficiency pairs
// (MCTrigEffPairSel::Step3PairSelection(true), RECO quantities) -- the population of the 4
// reference tables cells_{pTH70_125,pTH125_300}_{m_lt_10,m_sig}.csv of
// mc_pthat_slice_mass_statistics.cxx (reproduced here and checked against them).
//
// NEW: truth-only test samples of the request (1 000 / 100 events per config), processed by the
// ntuple processing in its external-truth raw-counts mode (PythiaTruthAnalysis,
// truth_input_override + truth_raw_counts). They have no reco, so the projection goes through a
// TRUTH-level rate in the same cell:
//
//   n_new(c) = n_ref,reco(c) x [N_req * nu_new(c)] / [N_ref * nu_ref(c)]
//   nu(c)    = truth pairs per processed event in cell c (truth pair pT, |truth eta^pair|, truth
//              mass window, truth sign tree), with the TRUTH analogue of every Step-3 leg cut:
//              truth fiducial (pT>4.5, |eta|<2.4), forward low-pT veto, gap windows on truth q*eta,
//              |truth eta^pair| < 2.2 -- the SAME helper strings, evaluated on truth columns.
//   boost(c) = [n_ref,reco(c) + n_new(c)] / n_ref,reco(c)       (per pT-hat slice)
//
// Assumption (stated, not tested here): the truth -> Step-3-reco efficiency of a cell is the same
// in the filtered and the unfiltered sample -- the cell is kinematic, and the filter acts on truth.
// Within ONE slice the old and new samples add as raw counts (in the filter-pass region a combined
// sample has one weight per slice). Counts are RAW and PER SLICE -- never summed across slices.
//
// Usage: root -l -b -q 'mc_request_stats_boost.cxx+("<dir of new-mode outputs>","<out dir>")'
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TSystem.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../../MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "../../../MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "../../../Utilities/MCTrigEffPairSelection.h"
#include "../../../Utilities/PairTrigEffEvaluator.h"
#include "dr_correction_sample_cfg.h"
#include "../../../Utilities/MCRequestLooseMassFilter.h"

namespace {

std::string F(const char* f, double a) { char b[256]; snprintf(b, sizeof b, f, a); return b; }

// The request: one entry per test configuration. kn = PythiaAlgCoreT kinRanges index (the per-kn
// tree index of the reference file). N_req = the event count requested from production.
struct NewCfg { std::string label; int kn; double n_req; std::string group; };
const std::vector<NewCfg> kNew = {
    {"5p36_SI_m10_pT45", 4, 1.0e6, "SI"}, {"5p36_SI_m10_pT45", 5, 1.0e6, "SI"},
    {"5p36_SS_m10",      2, 0.5e6, "SS"}, {"5p36_SS_m10_pT20", 3, 0.5e6, "SS"},
    {"5p36_SS_m10_pT20", 4, 0.5e6, "SS"}, {"5p36_SS_m10_pT20", 5, 0.5e6, "SS"},
};
const std::map<int, std::string> kSliceToken = {{2, "pTH24_40"}, {3, "pTH40_70"},
                                                {4, "pTH70_125"}, {5, "pTH125_300"}};

struct Meta { Long64_t nproc = 0, nbeam = 0; };
Meta ReadMeta(const std::string& file, int kn)
{
    std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));
    if (!f || f->IsZombie()) throw std::runtime_error("mc_request_stats_boost: cannot open " + file);
    TTree* t = dynamic_cast<TTree*>(f->Get("meta_tree_out"));
    const std::string bp = "nproc_kin" + std::to_string(kn) + "_beam0";
    const std::string bb = "nbeam_kin" + std::to_string(kn) + "_beam0";
    if (!t || t->GetEntries() < 1 || !t->GetBranch(bp.c_str()) || !t->GetBranch(bb.c_str()))
        throw std::runtime_error("mc_request_stats_boost: no usable meta_tree_out (" + bp + ") in "
                                 + file + " -- refusing: the per-event rate needs N_proc.");
    Meta m;
    t->SetBranchAddress(bp.c_str(), &m.nproc);
    t->SetBranchAddress(bb.c_str(), &m.nbeam);
    t->GetEntry(0);
    if (m.nproc <= 0) throw std::runtime_error("mc_request_stats_boost: " + bp + " <= 0 in " + file);
    return m;
}

// Aliases so the Step-3 helper strings (written on m1_pt, m1_eta, m1_charge, m1_truth_*, pair_eta)
// evaluate on TRUTH columns.
ROOT::RDF::RNode TruthAliases(ROOT::RDF::RNode d)
{
    return d.Alias("m1_pt", "m1.truth_pt").Alias("m1_eta", "m1.truth_eta")
            .Alias("m1_charge", "m1.truth_charge")
            .Alias("m1_truth_pt", "m1.truth_pt").Alias("m1_truth_eta", "m1.truth_eta")
            .Alias("m2_pt", "m2.truth_pt").Alias("m2_eta", "m2.truth_eta")
            .Alias("m2_charge", "m2.truth_charge")
            .Alias("m2_truth_pt", "m2.truth_pt").Alias("m2_truth_eta", "m2.truth_eta")
            .Alias("t_pair_eta", "truth_pair_eta")
            .Alias("c_pair_pt", "truth_pair_pt").Alias("c_pair_eta", "truth_pair_eta")
            .Alias("c_minv", "truth_minv");
}
std::string TruthStep3Analogue()
{
    // FiducialGapCut() hard-wires the column `pair_eta`, which already exists (RECO) in these
    // trees and cannot be aliased; so its three terms are rebuilt from the SAME ParamsSet
    // primitives it uses, with the pair term on the truth column. Mirror of
    // MCTrigEffPairSel::FiducialGapCut() -- keep in lockstep.
    return MCTrigEffPairSel::TruthFiducial() + " && " + MCTrigEffPairSel::ForwardLowPtVeto()
         + " && " + ParamsSet::FiducialGapCutExpr("m1_charge * m1_eta")
         + " && " + ParamsSet::FiducialGapCutExpr("m2_charge * m2_eta")
         + " && " + ParamsSet::PairFiducialEtaCutExpr("t_pair_eta");
}

struct Range { std::string token; double lo, hi; };
std::vector<Range> Ranges()
{
    const PairTrigEff::MassWindow& sig = PairTrigEff::Window("sig");
    return {{"m_lt_10", 0.0, MCRequest::kLooseMassMax}, {"m_sig", sig.lo, sig.hi}};
}

// cell counts (pair pT x |eta^pair|) of one tree for one mass range
ROOT::RDF::RResultPtr<TH2D> BookCells(ROOT::RDF::RNode d, const Range& R, const std::string& name,
                                      const std::vector<double>& pe, const std::vector<double>& ee)
{
    return d.Filter(F("c_minv > %.10g", R.lo) + F(" && c_minv < %.10g", R.hi))
            .Define("abs_pair_eta", "fabs(c_pair_eta)")
            .Histo2D({name.c_str(), "", int(pe.size()) - 1, pe.data(), int(ee.size()) - 1, ee.data()},
                     "c_pair_pt", "abs_pair_eta");
}
} // namespace

void mc_request_stats_boost(const char* new_dir, const char* out_dir, bool use_tight_wp = true)
{
    PairTrigEff::CheckSignalWindowMirror();
    const DrCorrSample cfg = GetDrCorrSample("pp_full", use_tight_wp);
    const std::string ref_file = cfg.sample_dir
        + "muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_full.root";
    const std::vector<double> pe = PairTrigEff::PairPtEdges();
    const std::vector<double> ee = PairTrigEff::AbsEtaGroups().edges;
    const int first = PairTrigEff::FirstDeliveredPtBin();          // 3 highest cells: first..npt
    const int npt = int(pe.size()) - 1, neta = int(ee.size()) - 1;
    const std::vector<Range> ranges = Ranges();
    const std::string wp_col = MCTrigEffPairSel::PairWpBranch(use_tight_wp);
    gSystem->mkdir(out_dir, kTRUE);

    std::cout << "reference: " << ref_file << "\n  reco selection : "
              << MCTrigEffPairSel::Step3PairSelection(true)
              << "\n  truth analogue : " << TruthStep3Analogue() << std::endl;

    // ---- book everything, then run
    struct Key { std::string what; int kn; int sign; std::string range; std::string label; };
    std::map<std::string, ROOT::RDF::RResultPtr<TH2D>> H;
    std::vector<std::unique_ptr<ROOT::RDataFrame>> keep;
    std::map<int, Meta> meta_ref;
    std::map<std::string, Meta> meta_new;   // label_kn
    auto hkey = [](const std::string& w, int kn, int s, const std::string& r, const std::string& l) {
        return w + "|" + std::to_string(kn) + "|" + std::to_string(s) + "|" + r + "|" + l; };

    for (int kn : {2, 3, 4, 5}) {
        meta_ref[kn] = ReadMeta(ref_file, kn);
        for (int s : {1, 2}) {
            const std::string tree = "muon_pair_tree_kin" + std::to_string(kn) + "_sign" + std::to_string(s);
            keep.emplace_back(std::make_unique<ROOT::RDataFrame>(tree, ref_file));
            MCTrigEffPairSel::RequirePairWpColumn(keep.back()->GetColumnNames(), wp_col, ref_file + ":" + tree);
            ROOT::RDF::RNode base = *keep.back();
            // reco Step-3 (exactly as mc_pthat_slice_mass_statistics.cxx)
            ROOT::RDF::RNode reco = base.Alias("pair_wp", wp_col)
                .Alias("m1_pt", "m1.pt").Alias("m1_eta", "m1.eta").Alias("m1_charge", "m1.charge")
                .Alias("m1_truth_pt", "m1.truth_pt").Alias("m1_truth_eta", "m1.truth_eta")
                .Alias("m2_pt", "m2.pt").Alias("m2_eta", "m2.eta").Alias("m2_charge", "m2.charge")
                .Alias("m2_truth_pt", "m2.truth_pt").Alias("m2_truth_eta", "m2.truth_eta")
                .Alias("c_pair_pt", "pair_pt").Alias("c_pair_eta", "pair_eta").Alias("c_minv", "minv")
                .Filter(MCTrigEffPairSel::Step3PairSelection(true));
            ROOT::RDF::RNode tru = TruthAliases(base).Filter(TruthStep3Analogue());
            for (const auto& R : ranges) {
                H[hkey("ref_reco", kn, s, R.token, "")]  = BookCells(reco, R, "rr" + tree + R.token, pe, ee);
                H[hkey("ref_truth", kn, s, R.token, "")] = BookCells(tru, R, "rt" + tree + R.token, pe, ee);
            }
        }
    }
    for (const auto& c : kNew) {
        const std::string f = std::string(new_dir) + "/muon_pairs_pythia_truth_" + c.label + "_kn"
                            + std::to_string(c.kn) + "_no_data_resonance_cuts.root";
        if (gSystem->AccessPathName(f.c_str())) { std::cerr << "MISSING " << f << " -- skipped" << std::endl; continue; }
        meta_new[c.label + "_" + std::to_string(c.kn)] = ReadMeta(f, c.kn);
        for (int s : {1, 2}) {
            const std::string tree = "muon_pair_tree_kin" + std::to_string(c.kn) + "_sign" + std::to_string(s);
            keep.emplace_back(std::make_unique<ROOT::RDataFrame>(tree, f));
            ROOT::RDF::RNode tru = TruthAliases(*keep.back()).Filter(TruthStep3Analogue());
            for (const auto& R : ranges)
                H[hkey("new_truth", c.kn, s, R.token, c.label)] =
                    BookCells(tru, R, "nt" + c.label + tree + R.token, pe, ee);
        }
    }

    // ---- tables
    const char* sign_name[3] = {"", "same sign", "opposite sign"};
    std::ofstream csv(std::string(out_dir) + "/stats_boost.csv");
    csv << "# Projected statistics boost of the additional MC request, per pT-hat slice (RAW counts, never summed across slices).\n"
        << "# reference = pp24 fullsim FULL Step-3 pairs (" << (use_tight_wp ? "Tight" : "Medium") << "), " << ref_file << "\n"
        << "# n_new = n_ref_reco * (N_req*nu_new)/(N_ref*nu_ref); nu = truth pairs/event in the cell (truth analogue of Step 3)\n"
        << "# rel_err_new = 1/sqrt(truth pairs in the cell in the TEST sample) (test-sample statistics only)\n"
        << "slice,sign,mass,pair_pt_lo,pair_pt_hi,abs_eta_lo,abs_eta_hi,n_ref_reco,ref_truth_pairs,N_ref,config,N_req,test_truth_pairs,N_test,n_new_projected,rel_err_new,boost\n";
    for (int kn : {2, 3, 4, 5}) {
        for (int s : {1, 2}) for (const auto& R : ranges) {
            TH2D* rr = H[hkey("ref_reco", kn, s, R.token, "")].GetPtr();
            TH2D* rt = H[hkey("ref_truth", kn, s, R.token, "")].GetPtr();
            std::vector<const NewCfg*> contrib;
            for (const auto& c : kNew) if (c.kn == kn && meta_new.count(c.label + "_" + std::to_string(kn))) contrib.push_back(&c);
            std::cout << "\n== " << kSliceToken.at(kn) << ", " << sign_name[s] << ", " << R.token
                      << "  (N_ref = " << meta_ref[kn].nproc << ")\n";
            std::cout << "  cell [pT^pair x |eta^pair|]          n_ref  |";
            for (auto* c : contrib) std::cout << " n_new(" << c->group << ")";
            std::cout << " | boost\n";
            for (int ix = first; ix <= npt; ++ix) for (int iy = 1; iy <= neta; ++iy) {
                const double nr = rr->GetBinContent(ix, iy), tr = rt->GetBinContent(ix, iy);
                double nsum = 0;
                std::cout << F("  [%6.1f,", pe[ix - 1]) << F("%6.1f) x ", pe[ix]) << F("[%.1f,", ee[iy - 1])
                          << F("%.1f)", ee[iy]) << F("  %8.0f |", nr);
                for (auto* c : contrib) {
                    const Meta& mn = meta_new[c->label + "_" + std::to_string(kn)];
                    const double tn = H[hkey("new_truth", kn, s, R.token, c->label)]->GetBinContent(ix, iy);
                    // tr == 0: no reference truth pair in the cell -> ratio undefined (reported as nan)
                    const double nn = tr > 0 ? nr * (c->n_req * tn / mn.nproc) / (meta_ref[kn].nproc * tr / meta_ref[kn].nproc) : NAN;
                    nsum += std::isnan(nn) ? 0 : nn;
                    if (tn > 0) std::cout << F("  %9.0f", nn) << F(" (±%2.0f%%)", 100. / std::sqrt(tn));
                    else {   // no test pair in the cell: 68% CL one-sided Poisson upper limit, 1.14 pairs
                        const double ul = tr > 0 ? nr * (c->n_req * 1.14 / mn.nproc) / tr : NAN;
                        std::cout << F("  <%8.0f", ul) << "  (0 test)";
                    }
                    csv << kSliceToken.at(kn) << "," << (s == 1 ? "SS" : "OS") << "," << R.token << ","
                        << pe[ix - 1] << "," << pe[ix] << "," << ee[iy - 1] << "," << ee[iy] << "," << nr << ","
                        << tr << "," << meta_ref[kn].nproc << "," << c->label << "," << c->n_req << "," << tn << ","
                        << mn.nproc << "," << nn << "," << (tn > 0 ? 1. / std::sqrt(tn) : NAN) << ","   // 0 test pairs: n_new = 0, upper limit in stdout
                        << (nr > 0 ? (nr + nn) / nr : NAN) << "\n";
                }
                std::cout << F(" | %7.1f", nr > 0 ? (nr + nsum) / nr : NAN) << "\n";
            }
        }
    }
    std::cout << "\nwritten " << out_dir << "/stats_boost.csv" << std::endl;
}
