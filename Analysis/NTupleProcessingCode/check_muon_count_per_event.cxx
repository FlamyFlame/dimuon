// Compare per-event muon multiplicity between two pbpb23 skim datasets.
// Usage: root -l -b -q 'check_muon_count_per_event.cxx'
// Tree: HeavyIonD3PD; muon count inferred from muon_pt vector size.

#include "ROOT/RDataFrame.hxx"
#include "TChain.h"
#include <cstdio>
#include <string>
#include <vector>
#include <glob.h>

std::vector<std::string> glob_files(const std::string& pattern) {
    glob_t g;
    std::vector<std::string> result;
    if (glob(pattern.c_str(), GLOB_TILDE, nullptr, &g) == 0)
        for (size_t i = 0; i < g.gl_pathc; ++i) result.emplace_back(g.gl_pathv[i]);
    globfree(&g);
    return result;
}

void print_table(const std::string& label, long long n1, long long n2, long long n3, long long ngt3) {
    long long total = n1 + n2 + n3 + ngt3;
    printf("\n=== %s  (total events: %lld) ===\n", label.c_str(), total);
    printf("  %-12s  %10s  %8s\n", "n_muon", "count", "fraction");
    printf("  %-12s  %10lld  %8.4f%%\n", "1",   n1,   100.0*n1  /total);
    printf("  %-12s  %10lld  %8.4f%%\n", "2",   n2,   100.0*n2  /total);
    printf("  %-12s  %10lld  %8.4f%%\n", "3",   n3,   100.0*n3  /total);
    printf("  %-12s  %10lld  %8.4f%%\n", ">3",  ngt3, 100.0*ngt3/total);
}

void run_dataset(const std::string& pattern, const std::string& label) {
    auto files = glob_files(pattern);
    if (files.empty()) { printf("WARNING: no files matched: %s\n", pattern.c_str()); return; }
    printf("Loading %zu file(s) for %s\n", files.size(), label.c_str());

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame rdf("HeavyIonD3PD", files);

    auto nmu = rdf.Define("n_muon", [](const std::vector<float>& pt){ return (int)pt.size(); },
                          {"muon_pt"});

    auto c1   = nmu.Filter("n_muon == 1").Count();
    auto c2   = nmu.Filter("n_muon == 2").Count();
    auto c3   = nmu.Filter("n_muon == 3").Count();
    auto cgt3 = nmu.Filter("n_muon >  3").Count();

    print_table(label, *c1, *c2, *c3, *cgt3);
}

void check_muon_count_per_event() {
    run_dataset("~/usatlasdata/dimuon_data/pbpb_2023/data_pbpb23_part*.root", "NEW skim (usatlasdata)");
    run_dataset("~/dcachearea/dimuon_data/pbpb_2023/data_pbpb23_part*.root",  "OLD skim (dcachearea)");
}
