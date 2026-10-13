// Output check of the run/LB exclusion (docs/tracking/ntuple_rerun_pp24_pbpb2526_20261013.md §3a):
// in a combined muon_pairs file, count OS+SS pairs with lb == 0 (lbn branch not filled) and pairs in an excluded
// run/LB (PbPbBadRuns/PbPbBadLBRanges, or PPBadLBRanges for pp).  Prints "EXCL nlb0 nbad maxlb nruns";
// exits 1 if nlb0 or nbad is non-zero.
#include <set>
#include <string>
#include "ROOT/RDataFrame.hxx"
#include "TSystem.h"
#include "../NTupleProcessingCode/PbPbEventSelConfig.h"
#include "../NTupleProcessingCode/PPEventSelConfig.h"

void check_runlb_exclusion(const char* file, bool is_pbpb, int year) {
    const auto pb_runs = PbPbBadRuns(year); const auto pb_lbs = PbPbBadLBRanges(year);
    const auto pp_lbs  = PPBadLBRanges(year);
    auto bad = [&](UInt_t r, UInt_t lb) {
        return is_pbpb ? PbPbIsBadRunLB(pb_runs, pb_lbs, (int)r, lb) : PPIsBadRunLB(pp_lbs, (int)r, lb);
    };
    long nlb0 = 0, nbad = 0; unsigned maxlb = 0; std::set<unsigned> runs;
    for (const char* tree : {"muon_pair_tree_sign1", "muon_pair_tree_sign2"}) {
        ROOT::RDataFrame df(tree, file);
        auto r  = df.Take<UInt_t>("run_number");
        auto lb = df.Take<UInt_t>("lb");
        for (size_t i = 0; i < r->size(); ++i) {
            if ((*lb)[i] == 0) ++nlb0;
            if (bad((*r)[i], (*lb)[i])) ++nbad;
            if ((*lb)[i] > maxlb) maxlb = (*lb)[i];
            runs.insert((*r)[i]);
        }
    }
    printf("EXCL %ld %ld %u %zu\n", nlb0, nbad, maxlb, runs.size());
    gSystem->Exit((nlb0 || nbad) ? 1 : 0);
}
