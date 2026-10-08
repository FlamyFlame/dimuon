// Provenance check for the pair-pT proxy study: <R_pair> per pT-hat slice (kin trees), both samples. Output: plots/truth_level_studies/pair_pt_proxy_b_hadron/R_mean_tables.txt
#include "/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/Utilities/SignalAcceptanceDefinition.h"
void diag2(){
 auto sel = SignalAcceptanceDefinition::NumExpr();
 const char* files[2]={"/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/muon_pairs_pythia_5p36TeV_no_data_resonance_cuts.root","/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_full.root"};
 const char* nm[2]={"central","fullsim"};
 for(int s=0;s<2;s++){
  printf("== %s (columns: slice, sum w, <R>, <R> in pair pT 9-12.8, n)\n",nm[s]);
  for(int k=0;k<6;k++){
   ROOT::RDataFrame d0(Form("muon_pair_tree_kin%d_sign2",k), files[s]);
   auto df=d0.Filter(sel).Define("b_pt",[](const std::vector<float>&v){return v.empty()?0.f:v[0]/1000.f;},{"m1_last_b_hadron_prt_pt_eta_phi_m"}).Filter("b_pt>0").Define("R","truth_pair_pt/b_pt").Define("wR","weight*R");
   auto lo=df.Filter("truth_pair_pt<12.8");
   auto sw=df.Sum("weight"); auto swr=df.Sum("wR"); auto n=df.Count();
   auto slw=lo.Sum("weight"); auto slwr=lo.Sum("wR");
   printf(" kn%d sumw=%.5g <R>=%.4f  <R>(9-12.8)=%.4f n=%llu\n",k,*sw,*swr/ *sw,*slwr/ *slw,*n);
  }
 }
}
