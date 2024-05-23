#include <string.h>
#include <vector>

const int nMCmodes = 2;
const int nSigns = 2;
const int nDphi = 2;
const int nBatches = 6;

std::string mcdir = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/";
std::vector<std::string> sub_dirs = {"bb_full_sample/", "cc_full_sample/"};

std::vector<std::vector<std::string>> fnames = 
  {{"muon_pairs_mc_truth_bb_1-5.root", "muon_pairs_mc_truth_bb_6-10.root", "muon_pairs_mc_truth_bb_11-15.root", "muon_pairs_mc_truth_bb_16-20.root", "muon_pairs_mc_truth_bb_21-25.root", "muon_pairs_mc_truth_bb_26-30.root"},
   {"muon_pairs_mc_truth_cc_1-5.root", "muon_pairs_mc_truth_cc_6-10.root", "muon_pairs_mc_truth_cc_11-15.root", "muon_pairs_mc_truth_cc_16-20.root", "muon_pairs_mc_truth_cc_21-25.root", "muon_pairs_mc_truth_cc_26-30.root"}};
std::vector<std::string> mcmodes = {"_bb","_cc"};
std::vector<std::string> signs = {"_sign1", "_sign2"};
std::vector<std::string> dphis = {"_near", "_away"};

void hist_helper(TH1* h, std::string title){
  h->SetStats(0);
  h->SetTitle(title.c_str());
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.6);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(38);
  h->GetYaxis()->SetLabelOffset(0.01);
  // h->GetYaxis()->SetTitleFont(43);
  // h->GetYaxis()->SetTitleSize(32);
  // h->GetYaxis()->SetTitleOffset(2);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(38);
  h->GetXaxis()->SetLabelOffset(0.01);
  // h->GetXaxis()->SetTitleFont(43);  
  // h->GetXaxis()->SetTitleSize(32);
  // h->GetXaxis()->SetTitleOffset(1);
}
