// Build data_pbpb23_part6.root = the events of PbPb23 runs 463315 / 462759 that the May-2026 skim LOST
// (JEDI --mergeOutput left finished job outputs unmerged; tracking doc data_per_run_signal_yield_qa.md Round 3/4).
//
// Input: the re-skim of the two runs (task 53029101, all events of both runs) and the May parts that hold them
// (462759 -> part2, 463315 -> part3). Output: re-skim events whose (RunNumber, eventNumber) is NOT in the May part,
// i.e. event-level dedup, independent of which job outputs were orphaned. The extra re-skim branch muon_match_L1MU3V
// is dropped so the schema equals part1-5 (common branches proven identical on overlapping events).
//
// Aborts (no output) on: re-skim internal duplicates, a re-skim run outside the two, or a May event missing from the
// re-skim (re-skim then incomplete, itself hit by the merge loss).
//
// root -l -b -q 'make_pbpb23_part6.C("rerun_dir/*.root", "out.root")'
#include <set>
#include <map>

static const char* kLGD = "root://dcgftp.usatlas.bnl.gov:1096//pnfs/usatlas.bnl.gov/LOCALGROUPDISK/rucio/user/yuhang/";

std::set<ULong64_t> may_events(const TString& path, UInt_t run) {
    TFile* f = TFile::Open(path);
    TTree* t = (TTree*)f->Get("HeavyIonD3PD");
    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("RunNumber", 1);
    t->SetBranchStatus("eventNumber", 1);
    UInt_t r; ULong64_t e;
    t->SetBranchAddress("RunNumber", &r);
    t->SetBranchAddress("eventNumber", &e);
    std::set<ULong64_t> s;
    Long64_t n = 0;
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        if (r != run) continue;
        ++n;
        if (!s.insert(e).second) { printf("FATAL: duplicate %u %llu in %s\n", r, e, path.Data()); gSystem->Exit(1); }
    }
    printf("May %s run %u: %lld events\n", path.Data(), run, n);
    delete f;
    return s;
}

void make_pbpb23_part6(const char* rerun_glob, const char* out_path,
                       const char* part2 = nullptr, const char* part3 = nullptr) {
    TString p2 = part2 ? TString(part2) : TString(kLGD) + "9a/b4/data_pbpb23_part2.root";
    TString p3 = part3 ? TString(part3) : TString(kLGD) + "be/82/data_pbpb23_part3.root";
    std::map<UInt_t, std::set<ULong64_t>> may;
    may[462759] = may_events(p2, 462759);
    may[463315] = may_events(p3, 463315);

    TChain c("HeavyIonD3PD");
    c.Add(rerun_glob);
    printf("re-skim: %d files, %lld entries\n", c.GetNtrees(), c.GetEntries());
    UInt_t r; ULong64_t e;
    c.SetBranchStatus("muon_match_L1MU3V", 0);
    c.SetBranchAddress("RunNumber", &r);
    c.SetBranchAddress("eventNumber", &e);

    // pass 1: decide which entries to keep, with every abort check, before any output is written
    std::map<UInt_t, std::set<ULong64_t>> seen;
    std::map<UInt_t, Long64_t> n_new, n_old;
    std::vector<Long64_t> keep;
    for (Long64_t i = 0; i < c.GetEntries(); ++i) {
        c.GetEntry(i);
        if (!may.count(r)) { printf("FATAL: re-skim run %u not in {462759, 463315}\n", r); gSystem->Exit(1); }
        if (!seen[r].insert(e).second) { printf("FATAL: re-skim duplicate %u %llu\n", r, e); gSystem->Exit(1); }
        if (may[r].count(e)) { ++n_old[r]; continue; }
        ++n_new[r];
        keep.push_back(i);
    }
    bool ok = true;
    for (auto& [run, s] : may) {
        Long64_t missing = 0;
        for (auto ev : s) if (!seen[run].count(ev)) ++missing;
        printf("run %u: re-skim %zu | May %zu | in both %lld | May-only %lld | NEW (-> part6) %lld\n",
               run, seen[run].size(), s.size(), n_old[run], missing, n_new[run]);
        if (missing) ok = false;
    }
    if (!ok) { printf("FATAL: May events missing from the re-skim -> re-skim incomplete, no output written\n"); gSystem->Exit(1); }

    // pass 2: write
    TFile out(out_path, "RECREATE");
    TTree* t = c.CloneTree(0);
    for (auto i : keep) { c.GetEntry(i); t->Fill(); }
    t->Write();
    printf("wrote %lld events, %d branches to %s\n", t->GetEntries(), t->GetListOfBranches()->GetEntries(), out_path);
    out.Close();
}
