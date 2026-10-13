# Filter-efficiency ESTIMATE for the new JOs, by emulating the new xAODDiLeptonMassFilter on the
# existing pp truth-only NTUPs (nPDF evgen e8599; every event already passed xAODMultiMuonFilter).
# eps_new ~= genFiltEff(AMI, PDF _pdf DSID of the same slice) * f,  f = fraction of NTUP events
# with >=1 truth mu pair: both legs pT>3.7 GeV |eta|<2.5, m<10 GeV, pT(mumu)>cut (SS if required).
import ROOT, sys
ROOT.gInterpreter.Declare("""
#include "ROOT/RVec.hxx"
bool passF(const ROOT::RVec<float>& pt1,const ROOT::RVec<float>& eta1,const ROOT::RVec<int>& ch1,
           const ROOT::RVec<float>& pt2,const ROOT::RVec<float>& eta2,const ROOT::RVec<int>& ch2,
           const ROOT::RVec<float>& m,const ROOT::RVec<float>& pt,double ptcut,bool ss,double sc){
  for(size_t i=0;i<m.size();++i){
    if(fabs(pt1[i])*sc<=3.7||fabs(pt2[i])*sc<=3.7||fabs(eta1[i])>2.5||fabs(eta2[i])>2.5) continue;
    if(ss && ch1[i]*ch2[i]<=0) continue;
    if(m[i]*sc>0 && m[i]*sc<10. && pt[i]*sc>ptcut) return true; }
  return false; }""")
d="/pnfs/usatlas.bnl.gov/users/yuhanguo/pythia_truth_full_sample/pythia_5p36TeV/"
ami="/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV/ami_info_PDF/ami_info_mc23_5p36TeV_Py8EG_A14_pp_hQCD_DiMu_pTH%s_pdf.txt"
cfgs=[("SI","70_125",45,False),("SI","125_300",45,False),("SS","24_40",20,True),("SS","40_70",20,True),("SS","70_125",20,True),("SS","125_300",20,True)]
maxev=int(sys.argv[1]) if len(sys.argv)>1 else 2000000
for g,sl,cut,ss in cfgs:
  eff=[float(l.split(":")[1]) for l in open(ami%sl) if l.startswith("genFiltEff ")][0]
  df=ROOT.RDataFrame("HeavyIonD3PD",d+"pythia_5p36TeV_pp_hQCD_DiMu_pTH%s.TRUTH0.NTUP.root"%sl).Range(maxev)
  # unit probe: muon pT in MeV if median > 1000
  sc=1e-3 if df.Range(2000).Mean("truth_muon_pt").GetValue()>100 else 1.
  n=df.Count(); k=df.Filter("passF(truth_mupair_pt1,truth_mupair_eta1,truth_mupair_ch1,truth_mupair_pt2,truth_mupair_eta2,truth_mupair_ch2,truth_mupair_m,truth_mupair_pt,%g,%s,%g)"%(cut,"true" if ss else "false",sc)).Count()
  f=k.GetValue()/n.GetValue()
  print("%s pTH%-8s genFiltEff=%.3e  f=%d/%d=%.4e  eps_est=%.3e  (unit scale %g)"%(g,sl,eff,k.GetValue(),n.GetValue(),f,eff*f,sc),flush=True)
