// check_event_duplicates.C -- prove that no (RunNumber, eventNumber) appears twice across a
// set of merged skim NTUPs (e.g. all parts of one year, including recovery parts).
//
//   root -l -b -q 'check_event_duplicates.C("/path/pbpb_2025/data_pbpb25_part*.root")'
//
// Strategy: a full hash set of ~3e8 events would need many GB, so work per RUN.  Pass 1
// reads only RunNumber from every file and records which files contain which runs.  Runs
// present in only ONE file cannot be duplicated across files (and within-file duplicates
// are impossible for a single skim job stream).  Runs present in >=2 files get an exact
// check in pass 2: read (RunNumber,eventNumber) from just those files, restricted to that
// run, into a std::unordered_set<ULong64_t>; any repeat is a duplicate and is printed.
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TObjArray.h>
#include <unordered_set>
#include <map>
#include <set>
#include <vector>
#include <cstdio>

static std::vector<TString> expand(const char* pattern){
   std::vector<TString> out;
   TString dir = gSystem->DirName(pattern), base = gSystem->BaseName(pattern);
   void* d = gSystem->OpenDirectory(dir);
   const char* e;
   TRegexp re(base, kTRUE);
   while(d && (e = gSystem->GetDirEntry(d))){ TString n(e); if(n.Index(re)==0 && n.EndsWith(".root")) out.push_back(dir+"/"+n); }
   gSystem->FreeDirectory(d);
   std::sort(out.begin(), out.end());
   return out;
}

void check_event_duplicates(const char* pattern, const char* tree="HeavyIonD3PD"){
   auto files = expand(pattern);
   printf("files (%zu):\n", files.size()); for(auto& f: files) printf("  %s\n", f.Data());

   // pass 1: run -> set of file indices, and entries per (file,run)
   std::map<UInt_t, std::set<size_t>> runFiles;
   std::map<std::pair<size_t,UInt_t>, Long64_t> nPerFileRun;
   Long64_t total=0;
   for(size_t i=0;i<files.size();++i){
      TFile* f=TFile::Open(files[i]); TTree* t=(TTree*)f->Get(tree);
      t->SetBranchStatus("*",0); t->SetBranchStatus("RunNumber",1);
      UInt_t run=0; t->SetBranchAddress("RunNumber",&run);
      Long64_t n=t->GetEntries(); total+=n;
      UInt_t last=0; Long64_t cnt=0;
      for(Long64_t k=0;k<n;++k){ t->GetEntry(k); runFiles[run].insert(i); nPerFileRun[{i,run}]++; }
      printf("  pass1 %s: %lld entries, %zu runs\n", gSystem->BaseName(files[i]), n, runFiles.size());
      f->Close();
   }
   printf("total entries %lld, distinct runs %zu\n", total, runFiles.size());

   // pass 2: exact check on runs shared by >=2 files
   Long64_t ndup=0, nchecked=0; int nshared=0;
   for(auto& kv: runFiles){
      if(kv.second.size()<2) continue;
      ++nshared;
      UInt_t run=kv.first;
      std::unordered_set<ULong64_t> seen; seen.reserve(4000000);
      std::map<ULong64_t,int> dupOwner;
      for(size_t i: kv.second){
         TFile* f=TFile::Open(files[i]); TTree* t=(TTree*)f->Get(tree);
         t->SetBranchStatus("*",0); t->SetBranchStatus("RunNumber",1); t->SetBranchStatus("eventNumber",1);
         UInt_t r=0; ULong64_t ev=0; t->SetBranchAddress("RunNumber",&r); t->SetBranchAddress("eventNumber",&ev);
         Long64_t n=t->GetEntries();
         for(Long64_t k=0;k<n;++k){ t->GetEntry(k); if(r!=run) continue; ++nchecked;
            if(!seen.insert(ev).second){ ++ndup; if(ndup<=10) printf("  DUPLICATE run %u event %llu (in %s)\n", run, ev, gSystem->BaseName(files[i])); }
         }
         f->Close();
      }
      printf("  run %u shared by %zu files: %zu unique events checked, %lld duplicates so far\n", run, kv.second.size(), seen.size(), ndup);
   }
   printf("\n=== RESULT: %d run(s) shared across files, %lld events checked exactly, %lld DUPLICATES ===\n", nshared, nchecked, ndup);
   printf(ndup==0 ? "=== NO DUPLICATE EVENTS ===\n" : "=== !!! DUPLICATES FOUND !!! ===\n");
}
