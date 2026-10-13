// Per-run signal-region OS pair yield from the NOMINAL muon-pair NTuples (no standalone re-derivation
// of the selection: the cuts are the ParamsSet expressions the crossx RDF fill applies, verbatim).
// Writes <outdir>/<tag>_per_run.csv : run, N_pairs, N_events, lumi (file unit), lumi_norm.
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>
#include <string>
#include "ROOT/RDataFrame.hxx"
#include "TSystem.h"
#include "MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "NTupleProcessingCode/PbPbEventSelConfig.h"

void per_run_signal_yield_table(const std::string& tag, const std::string& file, bool is_pbpb,
                                const std::string& outdir) {
    ROOT::EnableImplicitMT();
    const std::string sig =
        ParamsSet::SignalMinvCutExpr("minv") + " && " + ParamsSet::SignalPairPtCutExpr("pair_pt") + " && "
        + ParamsSet::FiducialGapCutExpr("m1.charge * m1.eta") + " && "
        + ParamsSet::FiducialGapCutExpr("m2.charge * m2.eta") + " && "
        + ParamsSet::PairFiducialEtaCutExpr("pair_eta");
    // sign2 = opposite sign (reference_sign_convention); WP = Tight (nominal); PbPb: centrality acceptance
    ROOT::RDataFrame df("muon_pair_tree_sign2", file);
    auto n_all = df.Count();
    ROOT::RDF::RNode n = df.Filter("pair_pass_tight", "tight");
    if (is_pbpb) {
        // data-quality exclusion (PbPbBadRuns / PbPbBadLBRanges): already applied at NTuple stage once the NTuple is
        // regenerated with it; repeated here so an older NTuple cannot leak excluded runs/LBs into the plot
        const int yr = std::stoi(tag.substr(tag.size() - 2));
        const auto bad_runs = PbPbBadRuns(yr); const auto bad_lbs = PbPbBadLBRanges(yr);
        n = n.Filter([bad_runs, bad_lbs](UInt_t r, UInt_t lb) { return !PbPbIsBadRunLB(bad_runs, bad_lbs, (int)r, lb); },
                     {"run_number", "lb"}, "good run/LB");
        const std::string ctr = "avg_centrality >= " + std::to_string(ParamsSet::ctrbins.front())
                              + " && avg_centrality < " + std::to_string(ParamsSet::ctrbins.back());
        n = n.Filter(ctr, "centrality acceptance");
    }
    n = n.Filter(sig, "signal");
    auto run = n.Take<UInt_t>("run_number");
    auto lb  = n.Take<UInt_t>("lb");
    auto bc  = n.Take<UInt_t>("bcid");
    auto ev  = n.Take<Int_t>("m1.ev_num");
    auto wgt = n.Take<double>("weight");
    auto rep = n.Report();
    std::cout << "[" << tag << "] pairs in tree: " << *n_all << "\n"; rep->Print();

    std::map<unsigned, long> npair; std::map<unsigned, std::set<std::tuple<unsigned,unsigned,int>>> evs;
    std::map<unsigned, double> sw;
    for (size_t i = 0; i < run->size(); ++i) {
        npair[(*run)[i]]++; evs[(*run)[i]].insert({(*lb)[i], (*bc)[i], (*ev)[i]}); sw[(*run)[i]] += (*wgt)[i];
    }
    std::ofstream o(outdir + "/" + tag + "_per_run.csv");
    o << "run,n_pairs,n_events,sum_weight\n";
    for (auto& kv : npair) o << kv.first << "," << kv.second << "," << evs[kv.first].size() << "," << sw[kv.first] << "\n";
}
