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
//   n_new(c) = N_req * nu_new(c) * eps(c),   eps(c) = sum_{kn=2..5, both signs} n_ref,reco(c)
//                                                   / sum_{kn=2..5, both signs} n_ref,truth(c)
//   (eps = truth -> Step-3-reco efficiency of the cell incl. migration, POOLED over slices and signs:
//   a per-slice, per-sign ratio is zero wherever the reference has no reco pair, i.e. exactly in the
//   cells the request exists to fill; tracking doc D5, revised 2026-09-25)
//   nu(c)    = truth pairs per processed event in cell c (truth pair pT, |truth eta^pair|, truth
//              mass window, truth sign tree), with the TRUTH analogue of every Step-3 leg cut:
//              truth fiducial (pT>4.5, |eta|<2.4), forward low-pT veto, gap windows on truth q*eta,
//              |truth eta^pair| < 2.2 -- the SAME helper strings, evaluated on truth columns.
//   boost(c) = [n_ref,reco(c) + n_new(c)] / n_ref,reco(c)       (per pT-hat slice, stdout only)
//
// Assumption (stated, not tested here): the truth -> Step-3-reco efficiency of a cell is the same
// in the filtered and the unfiltered sample -- the cell is kinematic, and the filter acts on truth.
// Within ONE slice the old and new samples add as raw counts (in the filter-pass region a combined
// sample has one weight per slice). stdout: per slice. OUTPUT (USER, 2026-10-02, D10): 4 tables
// projected_counts_{m_lt_10,m_sig}_{SS,OS}.csv, 3 pair-pT x 3 |eta^pair| cells, each cell the
// projected RAW pair count n_ref + n_new summed over the request's slices kn2-5 (a statistics
// count, not a yield: slices carry different sigma weights).
//
// DERIVED config (USER, 2026-10-04, D11): SS pTH125_300 m<10 pT^mumu>50 (0.2 M) has no truth test sample.
// It is the pT^mumu>20 SS sample restricted to the events whose FILTER-level SS pair (legs pT>3.7,
// |eta|<2.5, 0<m<10) has pT^mumu>50 -- flagged from the NTUP (the filter's own definition), matched to
// the processing output by m1.ev_num (validated pair by pair). nu = pairs in those events / N_subset.
// |eta^pair| dependence (user): in the TOP pair-pT bin (SS pT>20 statistics too low) the SS bin total
// is split by the OPPOSITE-SIGN truth-pair |eta^pair| fractions of the slice's test samples; in the
// other bins the SS cells are used as they are.
//
// Usage: root -l -b -q 'mc_request_stats_boost.cxx+("<dir of new-mode outputs>","<out dir>")'
#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TSystem.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
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

// 68.27% CL one-sided Poisson upper limit on the mean for 0 observed: -ln(1 - 0.6827)
const double kPoissonUL68_n0 = -std::log(1. - 0.6827);

std::string F(const char* f, double a) { char b[256]; snprintf(b, sizeof b, f, a); return b; }

// The request: one entry per test configuration. kn = PythiaAlgCoreT kinRanges index (the per-kn
// tree index of the reference file). N_req = the event count requested from production.
// Derived configs: src = label of the processed test sample, ntup_dir = its NTUP dir (under ntup_base),
// subset_pt = filter-level SS pair-pT threshold [GeV] defining the event subset (D11).
struct NewCfg { std::string label; int kn; double n_req; std::string group;
                std::string src = ""; std::string ntup_dir = ""; double subset_pt = -1.;
                std::string File() const { return src.empty() ? label : src; }
                bool Derived() const { return subset_pt > 0; } };
const std::vector<NewCfg> kNew = {
    {"5p36_SI_m10_pT45", 4, 1.0e6, "SI"}, {"5p36_SI_m10_pT45", 5, 1.0e6, "SI"},
    {"5p36_SS_m10",      2, 0.5e6, "SS"}, {"5p36_SS_m10_pT20", 3, 0.5e6, "SS"},
    {"5p36_SS_m10_pT20", 4, 0.5e6, "SS"}, {"5p36_SS_m10_pT20", 5, 0.5e6, "SS"},
    {"5p36_SS_m10_pT50", 5, 0.2e6, "SS50", "5p36_SS_m10_pT20", "run_same_sign_pT125_300", 50.},
};
// filter definition of the JO control file (xAODDiLeptonMassFilter on TruthMuons): legs pT > 3.7 GeV,
// |eta| < 2.5, 0 < m < kLooseMassMax, same sign
const double kFiltMuPt = 3.7, kFiltMuEta = 2.5;
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
    if (m.nproc != m.nbeam)   // a truncated processing would understate every per-event rate
        throw std::runtime_error("mc_request_stats_boost: " + file + " processed " + std::to_string(m.nproc)
                                 + " of " + std::to_string(m.nbeam) + " events for kn" + std::to_string(kn));
    return m;
}

// Per NTUP event (sorted file list = the processing input order): does it hold a filter-level SS pair
// with pT^mumu > pt_min? Also returns every pair pT [GeV] per event, for the ev_num match check.
struct Subset { std::vector<char> pass; std::vector<std::vector<float>> pair_pt; };
Subset FilterSubset(const std::string& dir, double pt_min)
{
    std::vector<std::string> files;
    FILE* p = gSystem->OpenPipe(("find " + dir + " -type f -name '*.root' ! -name '*log*' | LC_ALL=C sort").c_str(), "r");
    char buf[4096];
    while (p && fgets(buf, sizeof buf, p)) { std::string f(buf); f.erase(f.find_last_not_of("\n") + 1); files.push_back(f); }
    if (p) gSystem->ClosePipe(p);
    if (files.empty()) throw std::runtime_error("mc_request_stats_boost: no NTUP under " + dir);
    TChain ch("HeavyIonD3PD");
    for (const auto& f : files) ch.Add(f.c_str());
    std::vector<float> *pt1 = nullptr, *eta1 = nullptr, *pt2 = nullptr, *eta2 = nullptr, *m = nullptr, *pt = nullptr;
    std::vector<int> *c1 = nullptr, *c2 = nullptr;
    ch.SetBranchAddress("truth_mupair_pt1", &pt1); ch.SetBranchAddress("truth_mupair_eta1", &eta1);
    ch.SetBranchAddress("truth_mupair_ch1", &c1);  ch.SetBranchAddress("truth_mupair_pt2", &pt2);
    ch.SetBranchAddress("truth_mupair_eta2", &eta2); ch.SetBranchAddress("truth_mupair_ch2", &c2);
    ch.SetBranchAddress("truth_mupair_m", &m);     ch.SetBranchAddress("truth_mupair_pt", &pt);
    Subset S;
    for (Long64_t i = 0; i < ch.GetEntries(); ++i) {
        ch.GetEntry(i);
        bool pass = false; std::vector<float> pts;
        for (size_t k = 0; k < m->size(); ++k) {   // NTUP truth_mupair_* in MeV
            pts.push_back(pt->at(k) * 1e-3f);
            // filter convention (xAODDiLeptonMassFilter.cxx): legs pT >= MinPt, |eta| <= MaxEta; 0 < m < MaxMass; pT^mumu > Min
            if (std::fabs(pt1->at(k)) * 1e-3 < kFiltMuPt || std::fabs(pt2->at(k)) * 1e-3 < kFiltMuPt ||
                std::fabs(eta1->at(k)) > kFiltMuEta || std::fabs(eta2->at(k)) > kFiltMuEta) continue;
            const double mm = m->at(k) * 1e-3;
            if (c1->at(k) * c2->at(k) > 0 && mm > 0 && mm < MCRequest::kLooseMassMax && pt->at(k) * 1e-3 > pt_min) pass = true;
        }
        S.pass.push_back(pass); S.pair_pt.push_back(pts);
    }
    return S;
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

void mc_request_stats_boost(const char* new_dir, const char* out_dir, bool use_tight_wp = true,
                            const char* ntup_base = "/usatlas/u/yuhanguo/usatlasdata/evgen_high_pT_additional_test")
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
    std::map<std::string, std::shared_ptr<Subset>> subsets;   // derived configs, by label
    for (const auto& c : kNew) {
        const std::string f = std::string(new_dir) + "/muon_pairs_pythia_truth_" + c.File() + "_kn"
                            + std::to_string(c.kn) + "_no_data_resonance_cuts.root";
        if (gSystem->AccessPathName(f.c_str()))   // a missing config would silently drop out of the combined boost
            throw std::runtime_error("mc_request_stats_boost: missing new-sample file " + f);
        Meta mt = ReadMeta(f, c.kn);
        std::shared_ptr<Subset> sub;
        if (c.Derived()) {
            sub = std::make_shared<Subset>(FilterSubset(std::string(ntup_base) + "/" + c.ntup_dir, c.subset_pt));
            if ((Long64_t)sub->pass.size() != mt.nproc)
                throw std::runtime_error("mc_request_stats_boost: NTUP entries " + std::to_string(sub->pass.size())
                                         + " != processed events " + std::to_string(mt.nproc) + " for " + c.label);
            Long64_t nsub = 0; for (char b : sub->pass) nsub += b;
            std::cout << c.label << ": " << nsub << " of " << mt.nproc << " test events hold a filter-level SS pair with pT^mumu > "
                      << c.subset_pt << " GeV" << std::endl;
            mt.nproc = mt.nbeam = nsub;   // the derived sample's events
            subsets[c.label] = sub;
        }
        meta_new[c.label + "_" + std::to_string(c.kn)] = mt;
        for (int s : {1, 2}) {
            const std::string tree = "muon_pair_tree_kin" + std::to_string(c.kn) + "_sign" + std::to_string(s);
            keep.emplace_back(std::make_unique<ROOT::RDataFrame>(tree, f));
            ROOT::RDF::RNode d = *keep.back();
            if (sub) {
                // ev_num -> NTUP entry: every processed pair must exist (same truth pT) in that NTUP event
                auto bad = d.Filter([sub](int ev, float ppt) {
                    if (ev < 0 || ev >= (int)sub->pair_pt.size()) return true;
                    for (float x : sub->pair_pt[ev]) if (std::fabs(x - ppt) < 1e-3f * std::max(1.f, ppt)) return false;
                    return true; }, {"m1.ev_num", "truth_pair_pt"}).Count();
                if (*bad != 0) throw std::runtime_error("mc_request_stats_boost: " + std::to_string(*bad)
                                     + " processed pairs not found in their NTUP event (ev_num order mismatch) in " + f + ":" + tree);
                d = d.Filter([sub](int ev) { return (bool)sub->pass[ev]; }, {"m1.ev_num"});
            }
            ROOT::RDF::RNode tru = TruthAliases(d).Filter(TruthStep3Analogue());
            for (const auto& R : ranges)
                H[hkey("new_truth", c.kn, s, R.token, c.label)] =
                    BookCells(tru, R, "nt" + c.label + tree + R.token, pe, ee);
        }
    }
    // |eta^pair| shape source for the derived configs' TOP pair-pT bin (D11): opposite-sign truth pairs of
    // the slice's non-derived test samples, pooled
    std::set<int> os_eta_kn;
    for (const auto& c : kNew) {
        if (!c.Derived() || !os_eta_kn.insert(c.kn).second) continue;
        for (const auto& o : kNew) {
            if (o.kn != c.kn || o.Derived()) continue;
            const std::string f = std::string(new_dir) + "/muon_pairs_pythia_truth_" + o.File() + "_kn"
                                + std::to_string(o.kn) + "_no_data_resonance_cuts.root";
            keep.emplace_back(std::make_unique<ROOT::RDataFrame>("muon_pair_tree_kin" + std::to_string(o.kn) + "_sign2", f));
            ROOT::RDF::RNode tru = TruthAliases(*keep.back()).Filter(TruthStep3Analogue());
            for (const auto& R : ranges)
                H[hkey("os_eta_part", c.kn, 2, R.token, o.label)] = BookCells(tru, R, "oe" + o.label + R.token, pe, ee);
        }
    }
    // OS |eta^pair| counts of the top pair-pT bin, pooled over the slice's test samples: [kn][range][iy]
    auto OsEta = [&](int kn, const std::string& r, int iy) {
        double n = 0;
        for (const auto& o : kNew)
            if (o.kn == kn && !o.Derived()) n += H[hkey("os_eta_part", kn, 2, r, o.label)]->GetBinContent(npt, iy);
        return n; };

    // ---- pooled truth -> Step-3-reco efficiency per (mass range, cell), all slices, both signs
    std::map<std::string, std::unique_ptr<TH2D>> eps_reco, eps_truth;
    for (const auto& R : ranges) {
        for (int kn : {2, 3, 4, 5}) for (int s : {1, 2}) {
            TH2D* rr = H[hkey("ref_reco", kn, s, R.token, "")].GetPtr();
            TH2D* rt = H[hkey("ref_truth", kn, s, R.token, "")].GetPtr();
            if (!eps_reco[R.token]) { eps_reco[R.token].reset((TH2D*)rr->Clone(("er" + R.token).c_str())); eps_reco[R.token]->SetDirectory(nullptr);
                                      eps_truth[R.token].reset((TH2D*)rt->Clone(("et" + R.token).c_str())); eps_truth[R.token]->SetDirectory(nullptr); }
            else { eps_reco[R.token]->Add(rr); eps_truth[R.token]->Add(rt); }
        }
    }
    auto Eps = [&](const std::string& r, int ix, int iy) {
        const double t = eps_truth[r]->GetBinContent(ix, iy);
        return t > 0 ? eps_reco[r]->GetBinContent(ix, iy) / t : NAN; };

    // ---- tables
    // stdout: per slice, per-config projection n_new (0 test pairs -> n_new = 0 AND a 68% CL upper
    // limit) and the combined boost over the slice's configs (an upper value where configs had 0 test
    // pairs). Files: 4 projected-count tables, below.
    const char* sign_name[3] = {"", "same sign", "opposite sign"};
    for (const auto& R : ranges)   // every cell must have a reference efficiency -- checked BEFORE any output is written
        for (int ix = first; ix <= npt; ++ix) for (int iy = 1; iy <= neta; ++iy)
            if (std::isnan(Eps(R.token, ix, iy)))
                throw std::runtime_error("mc_request_stats_boost: eps undefined (no reference truth pair) in "
                                         + R.token + " cell " + std::to_string(ix) + "," + std::to_string(iy));
    // D11 checks. (a) An SS Step-3 pair in a delivered cell (pT >= 52 GeV, legs 4.5/2.4, m<10) is itself
    // a filter-level SS pair above the subset threshold, so the subset must hold ALL such pairs of the
    // source sample. (b) Print the OS |eta^pair| shape used in the top bin.
    for (const auto& c : kNew) {
        if (!c.Derived()) continue;
        for (const auto& R : ranges) {
            TH2D* hs = H[hkey("new_truth", c.kn, 1, R.token, c.label)].GetPtr();
            TH2D* hf = H[hkey("new_truth", c.kn, 1, R.token, c.File())].GetPtr();
            for (int ix = first; ix <= npt; ++ix) for (int iy = 1; iy <= neta; ++iy)
                if (pe[ix - 1] >= c.subset_pt && hs->GetBinContent(ix, iy) != hf->GetBinContent(ix, iy))
                    throw std::runtime_error("mc_request_stats_boost: subset lost SS pairs in a cell above the threshold ("
                                             + c.label + ", " + R.token + ")");
            std::cout << c.label << " " << R.token << ": SS test pairs per top-3 bin (subset)";
            for (int ix = first; ix <= npt; ++ix) { double t = 0; for (int iy = 1; iy <= neta; ++iy) t += hs->GetBinContent(ix, iy); std::cout << " " << t; }
            std::cout << "; top-bin OS |eta^pair| counts";
            for (int iy = 1; iy <= neta; ++iy) std::cout << " " << OsEta(c.kn, R.token, iy);
            std::cout << std::endl;
        }
    }
    // per (sign, mass range, cell): projected TOTAL raw pair count over the request's slices kn2-5
    // (USER, 2026-10-02, D10): n_ref,reco + sum n_new; 0-test configs enter with n_new = 0 (their
    // 68% CL upper values are listed in the CSV header).
    struct Tot { double ref = 0, nnew = 0, nnew_ul = 0; std::vector<std::string> zero; };
    std::map<std::string, Tot> tot;
    auto tkey = [](int s, const std::string& r, int ix, int iy) {
        return std::to_string(s) + "|" + r + "|" + std::to_string(ix) + "|" + std::to_string(iy); };
    for (int kn : {2, 3, 4, 5}) {
        std::vector<const NewCfg*> contrib;
        for (const auto& c : kNew) if (c.kn == kn) contrib.push_back(&c);
        for (int s : {1, 2}) for (const auto& R : ranges) {
            TH2D* rr = H[hkey("ref_reco", kn, s, R.token, "")].GetPtr();
            std::cout << "\n== " << kSliceToken.at(kn) << ", " << sign_name[s] << ", " << R.token
                      << "  (N_ref = " << meta_ref[kn].nproc << ")\n"
                      << "  cell [pT^pair x |eta^pair|]          n_ref  eps_pool |";
            for (auto* c : contrib) std::cout << "   n_new(" << c->group << ")        ";
            std::cout << "| boost (combined)\n";
            for (int ix = first; ix <= npt; ++ix) for (int iy = 1; iy <= neta; ++iy) {
                const double nr = rr->GetBinContent(ix, iy);
                const double eps = Eps(R.token, ix, iy);
                std::cout << F("  [%6.1f,", pe[ix - 1]) << F("%6.1f) x ", pe[ix]) << F("[%.1f,", ee[iy - 1])
                          << F("%.1f)", ee[iy]) << F("  %8.0f", nr) << F("  eps=%.3f |", eps);
                double nsum = 0, nsum_ul = 0; int nzero = 0;
                for (auto* c : contrib) {
                    const Meta& mn = meta_new.at(c->label + "_" + std::to_string(kn));
                    TH2D* hn = H[hkey("new_truth", kn, s, R.token, c->label)].GetPtr();
                    double tn = hn->GetBinContent(ix, iy), shape = 1.;   // shape: |eta| fraction applied to a bin total
                    double shape_ul = 0;   // > 0: the OS shape bin is empty -> UL of the cell from the shape statistics
                    if (c->Derived() && s == 1 && ix == npt) {
                        // D11: SS top-bin total, split by the OS |eta^pair| fractions of the slice's test samples
                        double bin_tot = 0, os_tot = 0;
                        for (int j = 1; j <= neta; ++j) { bin_tot += hn->GetBinContent(ix, j); os_tot += OsEta(kn, R.token, j); }
                        if (os_tot <= 0) throw std::runtime_error("mc_request_stats_boost: no OS test pair in the top bin for the |eta| shape ("
                                                                  + R.token + ", kn" + std::to_string(kn) + ")");
                        shape = OsEta(kn, R.token, iy) / os_tot;
                        if (shape == 0) shape_ul = kPoissonUL68_n0 / os_tot;
                        tn = bin_tot;   // statistics = the SS bin total (rel. error below: SS total only)
                    }
                    const bool zero_test = tn == 0 || shape_ul > 0;
                    const double nn = c->n_req * (tn / mn.nproc) * shape * eps;
                    const double ul = shape_ul > 0 ? c->n_req * (tn / mn.nproc) * shape_ul * eps
                                                   : c->n_req * (kPoissonUL68_n0 / mn.nproc) * shape * eps;
                    nsum += nn; nsum_ul += zero_test ? ul : nn; nzero += zero_test;
                    if (!zero_test) std::cout << F("  %9.0f", nn) << F(" (±%3.0f%%)", 100. / std::sqrt(tn));
                    else            std::cout << F("  <%8.0f", ul) << " (0 test)";
                    if (zero_test) tot[tkey(s, R.token, ix, iy)].zero.push_back(kSliceToken.at(kn) + " " + c->group
                                                                                + (shape_ul > 0 ? " (0 OS |eta| shape pairs)" : ""));
                }
                const double b = nr > 0 ? (nr + nsum) / nr : NAN, bul = nr > 0 ? (nr + nsum_ul) / nr : NAN;
                const bool all_zero = nzero == (int)contrib.size();
                if (nr == 0 && all_zero)              std::cout << F(" | n_ref = 0; n_new < %.0f (0 test)", nsum_ul);
                else if (nr == 0)                     std::cout << F(" | n_ref = 0; n_new = %.0f", nsum) << F(" (<= %.0f)", nsum_ul);
                else if (all_zero)                    std::cout << F(" | < %.1f (0 test)", bul);
                else if (nzero > 0)                   std::cout << F(" | %.1f", b) << F(" (<= %.1f)", bul);
                else                                  std::cout << F(" | %.1f", b);
                std::cout << "\n";
                Tot& T = tot[tkey(s, R.token, ix, iy)];
                T.ref += nr; T.nnew += nsum; T.nnew_ul += nsum_ul;
            }
        }
    }

    // ---- the 4 tables: (mass range) x (sign), 3 pair-pT rows x 3 |eta^pair| columns
    for (const auto& R : ranges) for (int s : {1, 2}) {
        const std::string path = std::string(out_dir) + "/projected_counts_" + R.token + "_" + (s == 1 ? "SS" : "OS") + ".csv";
        std::ofstream csv(path);
        if (!csv) throw std::runtime_error("mc_request_stats_boost: cannot open " + path + " for writing");
        csv << "# Additional pp24 Pythia MC request (ATLHI-666): PROJECTED pair count after the request, " << sign_name[s]
            << ", mass " << (R.lo > 0 ? F("%.2f < m < ", R.lo) : std::string("m < ")) << F("%.2f GeV", R.hi) << ".\n"
            << "# cell = n_ref + n_new: RAW Step-3 (" << (use_tight_wp ? "Tight" : "Medium") << ") pair count summed over pT-hat slices";
        for (int kn : {2, 3, 4, 5}) csv << " " << kSliceToken.at(kn);
        csv << "\n#   n_ref = existing pp24 fullsim FULL sample (" << ref_file << ")\n"
            << "#   n_new = N_req * nu_new * eps_pooled per requested config; nu_new = truth pairs/event of the grid test sample\n"
            << "#   in the cell (truth analogue of Step 3); eps_pooled = reference truth->Step-3-reco eff. Configs (N_req, N_test):\n#  ";
        for (const auto& c : kNew)
            csv << " {" << kSliceToken.at(c.kn) << " " << c.group << F(": %.0f, ", c.n_req)
                << meta_new.at(c.label + "_" + std::to_string(c.kn)).nproc << "}";
        csv << "\n"
            << "#   pTH125_300 SS50 = SS m<10 pT^mumu>50 (0.2 M), DERIVED from the pT>20 SS test sample: its events whose filter-level SS\n"
            << "#   pair has pT^mumu>50 (N_test above = that subset); in the top pair-pT bin its SS bin total is split by the OS |eta^pair|\n"
            << "#   fractions of the slice's test samples (D11).\n"
            << "#   (tracking doc additional_mc_request_jo_test.md D5 REVISED, D10, D11). Raw counts: slices carry different sigma weights.\n"
            << "# A config with 0 test pairs in a cell enters with n_new = 0 -> that cell is a LOWER value; listed below with the\n"
            << "#   68% CL upper value of the cell.\n";
        if (s == 2)
            csv << "# OS cells: the SS-filtered configs' contribution is approximate (an event with an OS pair here passes their filter\n"
                << "#   only if it also holds an SS pair); measured 2026-09-29 (before SS50): 0.9-2.4% of the SI contribution per table, up to 12.5% in single cells;\n"
                << "#   SS50 (D11) adds a further ~3-4% of the SI contribution in the m<10 cells it reaches.\n";
        for (int ix = first; ix <= npt; ++ix) for (int iy = 1; iy <= neta; ++iy) {
            const Tot& T = tot.at(tkey(s, R.token, ix, iy));
            if (T.zero.empty()) continue;
            csv << F("#   [%.1f,", pe[ix - 1]) << F("%.1f) x ", pe[ix]) << F("[%.1f,", ee[iy - 1]) << F("%.1f): 0 test pairs in", ee[iy]);
            for (const auto& z : T.zero) csv << " {" << z << "}";
            csv << F("; upper value %.0f\n", T.ref + T.nnew_ul);
        }
        csv << "# bin edges as used: pair pT [GeV]";
        for (int ix = first - 1; ix <= npt; ++ix) csv << F(" %.4f", pe[ix]);
        csv << "; |eta^pair|";
        for (double e : ee) csv << F(" %.4f", e);
        csv << "\npair_pT_GeV\\|pair_eta|";
        for (int iy = 1; iy <= neta; ++iy) csv << F(",%.1f-", ee[iy - 1]) << F("%.1f", ee[iy]);
        csv << "\n";
        for (int ix = first; ix <= npt; ++ix) {
            csv << F("%.1f-", pe[ix - 1]) << F("%.1f", pe[ix]);
            for (int iy = 1; iy <= neta; ++iy) { const Tot& T = tot.at(tkey(s, R.token, ix, iy)); csv << F(",%.0f", T.ref + T.nnew); }
            csv << "\n";
        }
        csv.close();
        if (!csv) throw std::runtime_error("mc_request_stats_boost: write to " + path + " FAILED (disk/quota?)");
        std::cout << "wrote " << path << std::endl;
    }
    std::cout << "\nNOTE: SS-filtered configs' contribution to OPPOSITE-SIGN cells is approximate (see CSV header).\n";
}
