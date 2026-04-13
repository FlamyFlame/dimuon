// trig_leg_check.C
//
// Produce per-leg trigger match counts and no-leg consistency tables
// for all dimuon triggers found in the input tree, at all dR cuts present.
//
// Usage:
//   root -l -b -q 'trig_leg_check.C("file.root")'
//   root -l -b -q 'trig_leg_check.C("f1.root,f2.root,...")'
//   root -l -b -q 'trig_leg_check.C("@filelist.txt")'          // one path per line
//   root -l -b -q 'trig_leg_check.C("file.root","MyTree")'     // custom tree name
//
// All input files must share the same directory.
// Output: <indir>/<stem>_trig_leg_check.md

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <iomanip>

void trig_leg_check(std::string input_arg, std::string tree_name = "HeavyIonD3PD") {

  // -----------------------------------------------------------------------
  // 1. Parse input files
  // -----------------------------------------------------------------------
  std::vector<std::string> files;
  if (!input_arg.empty() && input_arg[0] == '@') {
    std::ifstream fl(input_arg.substr(1));
    if (!fl) { std::cerr << "Cannot open file list: " << input_arg.substr(1) << "\n"; return; }
    std::string ln;
    while (std::getline(fl, ln)) {
      while (!ln.empty() && (ln.back()=='\r'||ln.back()==' ')) ln.pop_back();
      if (!ln.empty()) files.push_back(ln);
    }
  } else {
    std::stringstream ss(input_arg);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
      while (!tok.empty() && tok.front()==' ') tok.erase(tok.begin());
      while (!tok.empty() && (tok.back()==' '||tok.back()=='\r')) tok.pop_back();
      if (!tok.empty()) files.push_back(tok);
    }
  }
  if (files.empty()) { std::cerr << "No input files specified.\n"; return; }

  // -----------------------------------------------------------------------
  // 2. Validate same directory; build output path
  // -----------------------------------------------------------------------
  auto get_dir = [](const std::string& p) {
    size_t pos = p.rfind('/');
    return (pos == std::string::npos) ? std::string(".") : p.substr(0, pos);
  };
  auto get_stem = [](const std::string& p) {
    size_t pos = p.rfind('/');
    std::string base = (pos == std::string::npos) ? p : p.substr(pos+1);
    size_t dot = base.rfind('.');
    return (dot == std::string::npos) ? base : base.substr(0, dot);
  };

  std::string outdir = get_dir(files[0]);
  for (size_t i = 1; i < files.size(); i++) {
    if (get_dir(files[i]) != outdir) {
      std::cerr << "Error: all input files must be in the same directory.\n"
                << "  " << files[0] << "\n  " << files[i] << "\n";
      return;
    }
  }

  std::string outstem;
  if (files.size() == 1) {
    outstem = get_stem(files[0]);
  } else {
    outstem = get_stem(files[0]);
    for (size_t i = 1; i < files.size(); i++) {
      std::string s = get_stem(files[i]);
      size_t k = 0;
      while (k < outstem.size() && k < s.size() && outstem[k] == s[k]) k++;
      outstem = outstem.substr(0, k);
    }
    // strip trailing underscores/dashes; fall back if prefix is too short
    while (!outstem.empty() && (outstem.back()=='_'||outstem.back()=='-')) outstem.pop_back();
    if (outstem.size() < 3)
      outstem = "trig_leg_check_" + std::to_string(files.size()) + "files";
  }
  std::string outpath = outdir + "/" + outstem + "_trig_leg_check.md";

  // -----------------------------------------------------------------------
  // 3. Discover (trigger, dR) combinations from first file
  // -----------------------------------------------------------------------
  TFile* f0 = TFile::Open(files[0].c_str());
  if (!f0 || f0->IsZombie()) { std::cerr << "Cannot open: " << files[0] << "\n"; return; }
  TTree* t0 = (TTree*)f0->Get(tree_name.c_str());
  if (!t0) { std::cerr << "Tree \"" << tree_name << "\" not found in " << files[0] << "\n"; return; }

  const std::string PRE = "dimuon_b_";
  const std::string TAG = "_mu1passLeg1_dR_";

  // ordered list of (trigger_full_name, dr_str)
  std::vector<std::pair<std::string,std::string>> trig_dr_list;
  {
    TObjArray* brs = t0->GetListOfBranches();
    for (int i = 0; i < brs->GetEntries(); i++) {
      std::string bname = brs->At(i)->GetName();
      if (bname.size() <= PRE.size()) continue;
      if (bname.substr(0, PRE.size()) != PRE) continue;
      auto tpos = bname.find(TAG);
      if (tpos == std::string::npos) continue;
      std::string trig   = bname.substr(PRE.size(), tpos - PRE.size());
      std::string dr_str = bname.substr(tpos + TAG.size());

      // verify all 8 per-leg branches and no-leg branch exist
      bool ok = true;
      for (const char* mu : {"mu1","mu2"}) {
        for (const char* leg : {"Leg1","Leg2"}) {
          std::string b = PRE + trig + "_" + mu + "pass" + leg + "_dR_" + dr_str;
          if (!t0->GetBranch(b.c_str())) { ok = false; break; }
        }
        if (!ok) break;
      }
      if (ok && !t0->GetBranch((PRE + trig + "_" + dr_str).c_str())) ok = false;
      if (ok) trig_dr_list.push_back({trig, dr_str});
    }
  }
  f0->Close();

  if (trig_dr_list.empty()) {
    std::cerr << "No complete per-leg dimuon trigger branch sets found in " << files[0] << "\n";
    return;
  }

  std::cout << "Found " << trig_dr_list.size() << " (trigger, dR) combinations:\n";
  for (auto& p : trig_dr_list) std::cout << "  " << p.first << "  dR_" << p.second << "\n";

  // -----------------------------------------------------------------------
  // 4. Build TChain and set branch addresses
  // -----------------------------------------------------------------------
  TChain chain(tree_name.c_str());
  for (auto& fp : files) chain.Add(fp.c_str());

  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("muon_pair_muon1_index", 1);

  std::vector<int>* idx1 = nullptr;
  chain.SetBranchAddress("muon_pair_muon1_index", &idx1);

  struct BranchSet {
    std::vector<bool> *mu1L1=nullptr, *mu1L2=nullptr,
                      *mu2L1=nullptr, *mu2L2=nullptr, *noleg=nullptr;
  };
  std::vector<BranchSet> bsets(trig_dr_list.size());

  for (size_t ti = 0; ti < trig_dr_list.size(); ti++) {
    const auto& trig   = trig_dr_list[ti].first;
    const auto& dr_str = trig_dr_list[ti].second;
    auto& b = bsets[ti];

    auto enable = [&](const std::string& name, std::vector<bool>*& ptr) {
      chain.SetBranchStatus(name.c_str(), 1);
      chain.SetBranchAddress(name.c_str(), &ptr);
    };
    enable(PRE + trig + "_mu1passLeg1_dR_" + dr_str, b.mu1L1);
    enable(PRE + trig + "_mu1passLeg2_dR_" + dr_str, b.mu1L2);
    enable(PRE + trig + "_mu2passLeg1_dR_" + dr_str, b.mu2L1);
    enable(PRE + trig + "_mu2passLeg2_dR_" + dr_str, b.mu2L2);
    enable(PRE + trig + "_" + dr_str,                b.noleg);
  }

  // -----------------------------------------------------------------------
  // 5. Count loop
  // -----------------------------------------------------------------------
  struct Counts {
    long leg[4] = {};   // 0:mu1L1, 1:mu1L2, 2:mu2L1, 3:mu2L2
    long cons[4] = {};  // 0:both pass, 1:both fail, 2:combo pass & noleg fail, 3:combo fail & noleg pass
    long total = 0;
  };
  std::vector<Counts> counts(trig_dr_list.size());

  long nev = chain.GetEntries();
  for (long ev = 0; ev < nev; ev++) {
    chain.GetEntry(ev);
    int np = (int)idx1->size();
    for (size_t ti = 0; ti < trig_dr_list.size(); ti++) {
      auto& b = bsets[ti];
      auto& c = counts[ti];
      c.total += np;
      for (int ip = 0; ip < np; ip++) {
        if (b.mu1L1->at(ip)) c.leg[0]++;
        if (b.mu1L2->at(ip)) c.leg[1]++;
        if (b.mu2L1->at(ip)) c.leg[2]++;
        if (b.mu2L2->at(ip)) c.leg[3]++;
        bool combo = (b.mu1L1->at(ip) && b.mu2L2->at(ip)) ||
                     (b.mu1L2->at(ip) && b.mu2L1->at(ip));
        bool noleg = b.noleg->at(ip);
        if ( combo &&  noleg) c.cons[0]++;
        if (!combo && !noleg) c.cons[1]++;
        if ( combo && !noleg) c.cons[2]++;
        if (!combo &&  noleg) c.cons[3]++;
      }
    }
  }

  // -----------------------------------------------------------------------
  // 6. Write markdown
  // -----------------------------------------------------------------------
  std::ofstream md(outpath);
  md << "# Dimuon trigger leg-match statistics\n\n";
  md << "**Input files:**\n";
  for (auto& fp : files) md << "- `" << fp << "`\n";
  md << "\n**Tree:** `" << tree_name << "`  \n";
  md << "**Total events:** " << nev << "\n\n";

  for (size_t ti = 0; ti < trig_dr_list.size(); ti++) {
    const auto& trig   = trig_dr_list[ti].first;
    const auto& dr_str = trig_dr_list[ti].second;
    auto& c = counts[ti];

    std::string dr_disp = dr_str;
    std::replace(dr_disp.begin(), dr_disp.end(), '_', '.');

    md << "---\n\n## `" << trig << "` &nbsp; dR < " << dr_disp << "\n\n";
    md << "Total pairs: " << c.total << "\n\n";

    // leg counts
    md << "### Leg-pass counts\n\n";
    md << "| Muon | Leg | Count | % of pairs |\n";
    md << "|------|-----|------:|-----------:|\n";
    const char* mulab[] = {"mu1","mu1","mu2","mu2"};
    const char* leglab[] = {"Leg1","Leg2","Leg1","Leg2"};
    for (int i = 0; i < 4; i++) {
      char pct[16]; snprintf(pct, sizeof(pct), "%.3f%%", 100.0*c.leg[i]/c.total);
      md << "| " << mulab[i] << " | " << leglab[i] << " | "
         << c.leg[i] << " | " << pct << " |\n";
    }
    md << "\n";

    // consistency
    md << "### Consistency: `pair_pass` vs `(mu1L1&&mu2L2)||(mu1L2&&mu2L1)`\n\n";
    md << "| Outcome | Count | % of pairs |\n";
    md << "|---------|------:|-----------:|\n";
    const char* clab[] = {
      "leg-combo pass & pair pass (agree)",
      "leg-combo fail & pair fail (agree)",
      "leg-combo pass & pair **FAIL** (disagree)",
      "leg-combo **FAIL** & pair pass (disagree)"
    };
    for (int i = 0; i < 4; i++) {
      char pct[16]; snprintf(pct, sizeof(pct), "%.3f%%", 100.0*c.cons[i]/c.total);
      md << "| " << clab[i] << " | " << c.cons[i] << " | " << pct << " |\n";
    }
    md << "\n";
  }

  md.close();
  std::cout << "Written: " << outpath << "\n";
}
