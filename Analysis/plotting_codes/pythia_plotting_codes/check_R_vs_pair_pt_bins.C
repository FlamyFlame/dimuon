// Provenance check for the pair-pT proxy study (docs/tracking/single_b_pair_pt_proxy_b_hadron.md): <R_pair> per coarse pair-pT bin, both samples. Output: plots/truth_level_studies/pair_pt_proxy_b_hadron/R_mean_tables.txt
#include "/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/Utilities/SignalAcceptanceDefinition.h"
void diag(){
 ParamsSet pms; auto E = pms.pair_pt_coarse_bins; int nb=E.size()-1;
 auto sel = SignalAcceptanceDefinition::NumExpr();
 const char* files[2]={"/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/muon_pairs_pythia_5p36TeV_no_data_resonance_cuts.root","/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_full.root"};
 const char* nm[2]={"central","fullsim"};
 for(int s=0;s<2;s++){
  ROOT::RDataFrame d0("muon_pair_tree_sign2", files[s]);
  auto df=d0.Filter(sel).Define("b_pt",[](const std::vector<float>&v){return v.empty()?0.f:v[0]/1000.f;},{"m1_last_b_hadron_prt_pt_eta_phi_m"}).Filter("b_pt>0").Define("R","truth_pair_pt/b_pt").Define("wR","weight*R");
  auto hw=df.Histo1D({"hw","",nb,E.data()},"truth_pair_pt","weight");
  auto hwR=df.Histo1D({"hwr","",nb,E.data()},"truth_pair_pt","wR");
  double tot=hw->Integral(); printf("== %s tot=%g\n",nm[s],tot);
  for(int i=1;i<=nb;i++) printf(" pT[%5.1f,%5.1f) frac=%.4f <R>=%.4f\n",E[i-1],E[i],hw->GetBinContent(i)/tot,hwR->GetBinContent(i)/hw->GetBinContent(i));
 }
}
