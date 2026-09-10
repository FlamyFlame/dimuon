// check_skim_output.C -- sanity-check a TrigRates skim NTUP.
//
//   root -l -b -q 'scripts/check_skim_output.C("new.root","reference.root")'
//   root -l -b -q 'scripts/check_skim_output.C("new.root")'            // no comparison
//
// 1) Compares the branch list of <new> against a <reference> NTUP from a year that is
//    already validated, and prints branches present in one but not the other.
// 2) For every branch, reports the fraction of entries in which it is "filled":
//      - scalar        : value != 0 (and not NaN)
//      - fixed array   : any element != 0
//      - vector<T>     : size() > 0
//    A branch that exists but is 0/empty in 100% of events is flagged FILL=0 -- that is
//    the failure mode this macro exists to catch (a branch that silently never gets set
//    because a container key or a tool is missing for the new year).
//
// Trigger, muon scale-factor, ZDC and event-shape branches are additionally summarised
// in named groups so a missing group is obvious at a glance.

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {

std::set<std::string> branch_names(TTree* t) {
   std::set<std::string> out;
   TObjArray* bl = t->GetListOfBranches();
   for (int i = 0; i < bl->GetEntries(); ++i)
      out.insert(((TBranch*)bl->At(i))->GetName());
   return out;
}

// Fraction of entries in which the branch is non-trivially filled.
// Uses TTree::Draw with an expression that is robust for scalars, arrays and vectors.
double fill_fraction(TTree* t, const std::string& name, Long64_t nmax) {
   TBranch* b = t->GetBranch(name.c_str());
   if (!b) return -1.0;
   TLeaf* leaf = (TLeaf*)b->GetListOfLeaves()->At(0);
   const bool is_vector = (b->GetClassName() && std::string(b->GetClassName()).find("vector") != std::string::npos);

   std::string expr;
   if (is_vector)                       expr = "Length$(" + name + ")>0";
   else if (leaf && leaf->GetLen() > 1) expr = "Sum$(" + name + "!=0)>0";
   else                                 expr = name + "!=0";

   const Long64_t n = t->Draw(("(" + expr + ")").c_str(), "", "goff", nmax);
   if (n <= 0) return -1.0;
   const Double_t* v = t->GetV1();
   Long64_t filled = 0;
   for (Long64_t i = 0; i < n; ++i) if (v[i] != 0) ++filled;
   return double(filled) / double(n);
}

void group_report(TTree* t, const char* title, const std::vector<std::string>& pats,
                  const std::set<std::string>& names, Long64_t nmax) {
   std::cout << "\n--- " << title << " ---\n";
   int shown = 0;
   for (const std::string& nm : names) {
      bool match = false;
      for (const std::string& p : pats) if (nm.find(p) != std::string::npos) { match = true; break; }
      if (!match) continue;
      const double f = fill_fraction(t, nm, nmax);
      std::cout << std::left << std::setw(46) << nm << "  fill=";
      if (f < 0) std::cout << "n/a";
      else       std::cout << std::fixed << std::setprecision(4) << f;
      if (f == 0.0) std::cout << "   <== ALWAYS EMPTY";
      std::cout << "\n";
      ++shown;
   }
   if (!shown) std::cout << "  (no branches matched)   <== GROUP MISSING\n";
}

} // namespace

void check_skim_output(const char* new_file, const char* ref_file = nullptr,
                       Long64_t nmax = 200000, const char* tree_name = "HeavyIonD3PD") {
   TFile* fn = TFile::Open(new_file, "READ");
   if (!fn || fn->IsZombie()) { std::cout << "CANNOT OPEN " << new_file << "\n"; return; }
   TTree* tn = (TTree*)fn->Get(tree_name);
   if (!tn) { std::cout << "NO TREE '" << tree_name << "' IN " << new_file << "\n"; return; }

   std::cout << "\n==================== " << new_file << " ====================\n";
   std::cout << "tree     : " << tree_name << "\n";
   std::cout << "entries  : " << tn->GetEntries() << "\n";
   std::cout << "branches : " << tn->GetListOfBranches()->GetEntries() << "\n";

   std::set<std::string> nb = branch_names(tn);

   if (ref_file && *ref_file) {
      TFile* fr = TFile::Open(ref_file, "READ");
      if (fr && !fr->IsZombie()) {
         TTree* tr = (TTree*)fr->Get(tree_name);
         if (tr) {
            std::set<std::string> rb = branch_names(tr);
            std::cout << "\n=== branch-list diff vs reference " << ref_file << " ===\n";
            std::cout << "reference entries : " << tr->GetEntries()
                      << "   branches : " << rb.size() << "\n";
            int miss = 0, extra = 0;
            for (const std::string& s : rb) if (!nb.count(s)) { std::cout << "  MISSING in new : " << s << "\n"; ++miss; }
            for (const std::string& s : nb) if (!rb.count(s)) { std::cout << "  EXTRA   in new : " << s << "\n"; ++extra; }
            std::cout << "  => " << miss << " missing, " << extra << " extra\n";
            if (!miss && !extra) std::cout << "  => BRANCH LISTS IDENTICAL\n";
         }
         fr->Close();
      } else {
         std::cout << "\n(reference " << ref_file << " could not be opened -- skipping diff)\n";
      }
   }

   const Long64_t n = std::min<Long64_t>(nmax, tn->GetEntries());
   std::cout << "\n=== fill fractions (first " << n << " entries) ===\n";

   group_report(tn, "event info / vertex",   {"RunNumber","lbn","bcid","eventNumber","IntPerXing","vtx","Vtx"}, nb, n);
   group_report(tn, "muon kinematics",       {"mu_pt","mu_eta","mu_phi","mu_charge","mu_quality","mu_author","mu_type"}, nb, n);
   group_report(tn, "muon SF / efficiency corr", {"_SF","eff_corr","eff_SF","quality"}, nb, n);
   group_report(tn, "trigger decisions",     {"trig","Trig","HLT","L1_","hlt"}, nb, n);
   group_report(tn, "ZDC",                   {"zdc_"}, nb, n);
   group_report(tn, "event shape / FCal",    {"FCal","cent","Cent","Qvec","psi"}, nb, n);
   group_report(tn, "tracks",                {"trk_"}, nb, n);

   std::cout << "\n=== ALWAYS-EMPTY branches (full list) ===\n";
   int nempty = 0;
   for (const std::string& s : nb) {
      const double f = fill_fraction(tn, s, n);
      if (f == 0.0) { std::cout << "  " << s << "\n"; ++nempty; }
   }
   std::cout << "  => " << nempty << " of " << nb.size() << " branches are always empty/zero\n";
   std::cout << "\n==================== done ====================\n";
   fn->Close();
}
