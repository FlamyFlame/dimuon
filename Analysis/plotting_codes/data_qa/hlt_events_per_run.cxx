// Abnormal-run investigation (docs/tracking/data_per_run_signal_yield_qa.md): number of events firing one HLT chain,
// per run, from the raw skim NTUP (event-level trigger decision, no offline cut).  Used for pp24, which has no event
// selection beyond the bad-LB exclusion: N_sig/L = (N_HLT_2mu4/L) x (N_sig/N_HLT_2mu4).
// The pp bad-LB ranges PPBadLBRanges(pp_year) (pp24: run 488456 LB 84-220) are dropped, as in the NTuple stage.
// Writes <outcsv>: run,n_hlt
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "ROOT/RDataFrame.hxx"
#include "TChain.h"
#include "../../NTupleProcessingCode/PPEventSelConfig.h"

void hlt_events_per_run(const std::vector<std::string>& files, const std::string& trig_branch, const std::string& outcsv,
                        int pp_year = 24) {
    ROOT::EnableImplicitMT();
    TChain ch("HeavyIonD3PD");
    for (auto& f : files) ch.Add(f.c_str());
    ROOT::RDataFrame df(ch);
    const auto bad_lbs = PPBadLBRanges(pp_year);
    auto runs = df.Filter([](bool t) { return t; }, {trig_branch})
                  .Filter([bad_lbs](UInt_t r, UInt_t lb) { return !PPIsBadRunLB(bad_lbs, (int)r, lb); }, {"RunNumber", "lbn"})
                  .Take<UInt_t>("RunNumber");
    std::map<unsigned, long> n;
    for (auto r : *runs) n[r]++;
    std::ofstream o(outcsv);
    o << "run,n_hlt\n";
    for (auto& kv : n) o << kv.first << "," << kv.second << "\n";
    std::cout << "wrote " << outcsv << " (" << n.size() << " runs)\n";
}
