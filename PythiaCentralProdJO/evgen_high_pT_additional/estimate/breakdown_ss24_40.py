# SS pTH24_40: which cut costs what. Fraction of NTUP events (all passed the >=2mu pT>3.7 |eta|<2.5
# filter, AMI genFiltEff) having >=1 pair with both legs pT>3.7, |eta|<2.5 and the listed requirements.
import ROOT
ROOT.gInterpreter.Declare("""
#include "ROOT/RVec.hxx"
bool anyp(const ROOT::RVec<float>& pt1,const ROOT::RVec<float>& eta1,const ROOT::RVec<int>& ch1,
          const ROOT::RVec<float>& pt2,const ROOT::RVec<float>& eta2,const ROOT::RVec<int>& ch2,
          const ROOT::RVec<float>& m,const ROOT::RVec<float>& pt,int ss,double mmax,double ptmin){
  for(size_t i=0;i<m.size();++i){
    if(fabs(pt1[i])<=3700||fabs(pt2[i])<=3700||fabs(eta1[i])>2.5||fabs(eta2[i])>2.5) continue;
    if(ss==1 && ch1[i]*ch2[i]<=0) continue;
    if(ss==-1 && ch1[i]*ch2[i]>=0) continue;
    if(m[i]<mmax*1e3 && pt[i]>ptmin*1e3) return true; }
  return false; }""")
f="/pnfs/usatlas.bnl.gov/users/yuhanguo/pythia_truth_full_sample/pythia_5p36TeV/pythia_5p36TeV_pp_hQCD_DiMu_pTH24_40.TRUTH0.NTUP.root"
eff=1.821213e-4
df=ROOT.RDataFrame("HeavyIonD3PD",f)
args="truth_mupair_pt1,truth_mupair_eta1,truth_mupair_ch1,truth_mupair_pt2,truth_mupair_eta2,truth_mupair_ch2,truth_mupair_m,truth_mupair_pt"
cases=[("any sign, no m, no pT (=existing filter)",0,1e9,-1),("OS, no m, no pT",-1,1e9,-1),("SS, no m, no pT",1,1e9,-1),
       ("SS, m<10, no pT",1,10,-1),("SS, no m, pT>20",1,1e9,20),("SS, m<10, pT>20 (JO)",1,10,20),
       ("SS, m<10, pT>10",1,10,10),("SS, m<10, pT>15",1,10,15),("OS, m<10, pT>20",-1,10,20)]
n=df.Count(); cs=[df.Filter("anyp(%s,%d,%g,%g)"%(args,s,mm,pt)).Count() for _,s,mm,pt in cases]
N=n.GetValue(); print("NTUP events:",N)
for (lab,_,_,_),c in zip(cases,cs):
    k=c.GetValue(); e=eff*k/N; print("%-42s f=%7d/%d=%.3e  eps=%.2e  CPU/accepted=%.3g s"%(lab,k,N,k/N,e,0.036/e if e>0 else float('inf')))
