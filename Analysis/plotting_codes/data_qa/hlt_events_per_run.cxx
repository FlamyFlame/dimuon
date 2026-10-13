// Abnormal-run investigation (docs/tracking/data_per_run_signal_yield_qa.md): number of events firing one HLT chain,
// per run, from the raw skim NTUP (event-level trigger decision, no offline cut).  Used for pp24, which has no event
// selection: N_sig/L = (N_HLT_2mu4/L) x (N_sig/N_HLT_2mu4).
// Writes <outcsv>: run,n_hlt
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "ROOT/RDataFrame.hxx"
#include "TChain.h"

void hlt_events_per_run(const std::vector<std::string>& files, const std::string& trig_branch, const std::string& outcsv) {
    ROOT::EnableImplicitMT();
    TChain ch("HeavyIonD3PD");
    for (auto& f : files) ch.Add(f.c_str());
    ROOT::RDataFrame df(ch);
    auto runs = df.Filter([](bool t) { return t; }, {trig_branch}).Take<UInt_t>("RunNumber");
    std::map<unsigned, long> n;
    for (auto r : *runs) n[r]++;
    std::ofstream o(outcsv);
    o << "run,n_hlt\n";
    for (auto& kv : n) o << kv.first << "," << kv.second << "\n";
    std::cout << "wrote " << outcsv << " (" << n.size() << " runs)\n";
}
