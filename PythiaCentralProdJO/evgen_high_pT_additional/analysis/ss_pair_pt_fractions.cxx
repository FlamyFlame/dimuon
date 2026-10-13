// Why so few SS pairs in the top pair-pT cells? Test sample SS pTH125_300 (1000 evts).
#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/MuonObjectsParamsAndHelpers/FullSimSampleType.h"
#include "/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/MuonObjectsParamsAndHelpers/ParamsSet.h"
#include "/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/Utilities/MCTrigEffPairSelection.h"
#include "/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/Utilities/PairTrigEffEvaluator.h"
void ss_pair_pt_fractions(const char* ntup_glob, const char* pairfile, int kn) {
  const std::vector<double> E = {20, 30, 40, 52.2273, 74.2385, 105.5262, 150, 1e9};
  auto bin = [&](double x){ for (size_t i=0;i+1<E.size();++i) if (x>=E[i]&&x<E[i+1]) return int(i); return -1; };
  // ---- A: filter level (stable truth mu pT>3.7, |eta|<2.5), per event: LEADING qualifying SS m<10 pair
  TChain ch("HeavyIonD3PD"); ch.Add(ntup_glob);
  std::vector<float> *pt1=0,*eta1=0,*pt2=0,*eta2=0,*m=0,*pt=0; std::vector<int> *c1=0,*c2=0;
  ch.SetBranchAddress("truth_mupair_pt1",&pt1); ch.SetBranchAddress("truth_mupair_eta1",&eta1); ch.SetBranchAddress("truth_mupair_ch1",&c1);
  ch.SetBranchAddress("truth_mupair_pt2",&pt2); ch.SetBranchAddress("truth_mupair_eta2",&eta2); ch.SetBranchAddress("truth_mupair_ch2",&c2);
  ch.SetBranchAddress("truth_mupair_m",&m); ch.SetBranchAddress("truth_mupair_pt",&pt);
  std::vector<int> nev(E.size()-1,0), npair(E.size()-1,0); long N=ch.GetEntries(), nq=0, nmult=0;
  for (long i=0;i<N;++i){ ch.GetEntry(i); double lead=-1; int k_ok=0;
    for (size_t k=0;k<m->size();++k){
      if (fabs(pt1->at(k))*1e-3<=3.7||fabs(pt2->at(k))*1e-3<=3.7||fabs(eta1->at(k))>=2.5||fabs(eta2->at(k))>=2.5) continue;
      if (c1->at(k)*c2->at(k)<=0) continue;
      double mm=m->at(k)*1e-3, pp=pt->at(k)*1e-3; if(!(mm>0&&mm<10&&pp>20)) continue;
      ++k_ok; int b=bin(pp); if(b>=0) npair[b]++; lead=std::max(lead,pp);}
    if (k_ok) ++nq; if (k_ok>1) ++nmult; int b=bin(lead); if(b>=0) nev[b]++; }
  std::cout<<"A) filter level: events "<<N<<", with >=1 qualifying SS pair "<<nq<<", with >1: "<<nmult<<"\n   pair-pT bin : events(leading pair)  pairs\n";
  for (size_t i=0;i+1<E.size();++i) printf("   [%6.1f,%6.1f): %5d (%5.1f%%)  %5d\n",E[i],std::min(E[i+1],999.),nev[i],100.*nev[i]/N,npair[i]);
  // ---- B: ntuple-processing output (truth pairs, SS tree), analysis definitions
  ROOT::RDataFrame d("muon_pair_tree_kin"+std::to_string(kn)+"_sign1", pairfile);
  auto t = d.Alias("m1_pt","m1.truth_pt").Alias("m1_eta","m1.truth_eta").Alias("m1_charge","m1.truth_charge")
    .Alias("m1_truth_pt","m1.truth_pt").Alias("m1_truth_eta","m1.truth_eta")
    .Alias("m2_pt","m2.truth_pt").Alias("m2_eta","m2.truth_eta").Alias("m2_charge","m2.truth_charge")
    .Alias("m2_truth_pt","m2.truth_pt").Alias("m2_truth_eta","m2.truth_eta").Alias("t_pair_eta","truth_pair_eta")
    .Filter("truth_minv>0 && truth_minv<10");
  const std::string fid = MCTrigEffPairSel::TruthFiducial();
  const std::string step3 = fid+" && "+MCTrigEffPairSel::ForwardLowPtVeto()+" && "+ParamsSet::FiducialGapCutExpr("m1_charge * m1_eta")
     +" && "+ParamsSet::FiducialGapCutExpr("m2_charge * m2_eta")+" && "+ParamsSet::PairFiducialEtaCutExpr("t_pair_eta");
  auto hall = t.Histo1D({"a","",int(E.size())-1,E.data()},"truth_pair_pt");
  auto hfid = t.Filter(fid).Histo1D({"b","",int(E.size())-1,E.data()},"truth_pair_pt");
  auto hs3  = t.Filter(step3).Histo1D({"c","",int(E.size())-1,E.data()},"truth_pair_pt");
  std::cout<<"B) processing output, SS m<10 truth pairs: all | +truth fiducial (pT>4.5,|eta|<2.4) | +full Step-3 analogue\n";
  for (int i=1;i<=hall->GetNbinsX();++i) printf("   [%6.1f,%6.1f): %5.0f  %5.0f  %5.0f\n",E[i-1],std::min(E[i],999.),hall->GetBinContent(i),hfid->GetBinContent(i),hs3->GetBinContent(i));
}
