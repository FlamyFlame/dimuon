using namespace ROOT::VecOps;
// Mirrors NTupleProcessingCode single-muon selection:
//   DimuonDataAlgCoreT::PassCuts_DataCore (data, requireTight=true)
//   PythiaFullSimExtras::PassMuonMediumCuts + pass_tight       (MC)
// quality bits: 1=combined, 8=Medium, 16=Tight, 32=IDCuts, 256=MuonCuts  (305 = 1|16|32|256)
// pt is SIGNED in the ntuple (charge = sign); processing takes fabs()/1000 -> GeV
// ParamsSet.h: deltaP_overP_thrsh=0.12, d0cut=2mm, z0cut=2mm; turn_on_track_charge=false
// dpop_mode: 0 = |dP/P|<=0.12 (data code, fabs)   1 = dP/P<=0.12 (MC code, no fabs)
RVec<int> passmu(const RVec<float>&pt, const RVec<float>&eta, const RVec<int>&q,
                 const RVec<float>&dpop, const RVec<float>&d0, const RVec<float>&z0,
                 int dpop_mode){
  RVec<int> ok(pt.size(),0);
  for(size_t i=0;i<pt.size();++i){
    if((q[i]&305)!=305) continue;                       // combined & tight & IDCuts & MuonCuts
    if(std::fabs(eta[i])>2.4f) continue;                // |eta| <= 2.4
    if(std::fabs(pt[i])/1000.f < 4.f) continue;         // pt >= 4 GeV
    float dp = dpop[i];
    if(dpop_mode==0){ if(std::fabs(dp)>0.12f) continue; } else { if(dp>0.12f) continue; }
    if(std::fabs(d0[i])>=2.f) continue;                 // |d0| < 2 mm
    float z0s = std::fabs(z0[i]*std::sin(2.0*std::atan(std::exp(-eta[i]))));
    if(z0s>=2.f) continue;                              // |z0 sin(theta)| < 2 mm
    ok[i]=1;
  }
  return ok;
}
float subleadpt_sel(const RVec<float>&pt, const RVec<int>&ok){
  RVec<float> s;
  for(size_t i=0;i<pt.size();++i) if(ok[i]) s.push_back(std::fabs(pt[i])/1000.f);
  if(s.size()<2) return -1.f;
  std::sort(s.begin(),s.end(),std::greater<float>());
  return s[1];
}
