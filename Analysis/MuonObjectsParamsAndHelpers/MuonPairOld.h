#ifndef MuonPair_h
#define MuonPair_h

#include "Muon.h"
#include <iostream>
#include <string.h>
#include "TLorentzVector.h"
#include <math.h> 
#include <assert.h> 


// #include <fstream>
// #include "vector"

class MuonPair{
public:
	Muon m1;
	Muon m2;
  double weight;
  float pair_dPoverP; 
  float pt_lead;
  float pair_pt;
  float pair_eta;
  float pair_phi;
  float pair_y;
  float dpt;
  float deta;
  float dphi;
  float ptavg;
  float etaavg;
  float phiavg;
  float dr;
  float minv;
  float asym; // asymmetry := (pT_lead - pT_sublead) / (pT_lead + pT_sublead)
  float acop; // acoplanarity := (pi - |Dphi|) / pi
  bool same_sign;
  int avg_centrality;

	MuonPair(float pt1=0,float pt2=0,float eta1=5,float eta2=5,float phi1=5,float phi2=5,float dPoverP1=5,float dPoverP2=5);
	~MuonPair(){}
  // void InitPair();
  void InitPair(float pt1,float pt2,float eta1,float eta2,float phi1,float phi2,float dPoverP1,float dPoverP2);
	virtual void Sort();
  virtual void PairValueCalc();
  virtual void Update();
  virtual void UpdateShort();
  // void Sort();
  // void PairValueCalc();
  // void Update();
  // void UpdateShort();
};

MuonPair::MuonPair(float pt1=0,float pt2=0,float eta1=5,float eta2=5,float phi1=5,float phi2=5,float dPoverP1=5,float dPoverP2=5){
	weight = 1.;
  InitPair(pt1,pt2,eta1,eta2,phi1,phi2,dPoverP1,dPoverP2);
  // Update();
}

void MuonPair::InitPair(float pt1,float pt2,float eta1,float eta2,float phi1,float phi2,float dPoverP1,float dPoverP2){
  //initiaze the kinem1tics to unreasonable values to help with sanity check
  m1.pt = pt1;
  m2.pt = pt2;
  m1.eta = eta1;
  m2.eta = eta2;
  m1.phi = phi1;
  m2.phi = phi2;
  m1.dP_overP = dPoverP1;
  m2.dP_overP = dPoverP2;
}

void MuonPair::Sort(){
  ///Sort pt, eta, phi by pt

  float temppt, tempeta, tempphi, tempd0, tempz0, tempdP_overP;
  int tempind, tempcharge, tempquality;
  if (m1.pt < m2.pt){
    temppt = m1.pt;
    tempphi = m1.phi;
    tempeta = m1.eta;
    tempd0 = m1.d0;
    tempz0 = m1.z0;
    tempind = m1.ind;
    tempcharge = m1.charge;
    tempquality = m1.quality;
    tempdP_overP = m1.dP_overP;
  
    m1.pt = m2.pt;
    m1.eta = m2.eta;
    m1.phi = m2.phi;
    m1.d0 = m2.d0;
    m1.z0 = m2.z0;
    m1.ind = m2.ind;
    m1.charge = m2.charge;
    m1.quality = m2.quality;
    m1.dP_overP = m2.dP_overP;
  
    m2.pt = temppt;
    m2.eta = tempeta;
    m2.phi = tempphi;
    m2.d0 = tempd0;
    m2.z0 = tempz0;
    m2.ind = tempind;
    m2.charge = tempcharge;
    m2.quality = tempquality;
    m2.dP_overP = tempdP_overP;
  }
  // std::cout << m1.ev_num << ", parent class Sorting method called" << std::endl;
}

void MuonPair::PairValueCalc(){
  double PI=acos(-1.0);

  pt_lead = m1.pt;
  pair_dPoverP = sqrt(m1.dP_overP*m1.dP_overP + m2.dP_overP*m2.dP_overP);
  dpt = m1.pt - m2.pt;
  dphi = m1.phi-m2.phi;
  dphi = atan2(sin(dphi),cos(dphi));//fold dphi to [-pi,pi]
  deta = m1.eta-m2.eta;
  dr   = sqrt(dphi*dphi + deta*deta);
  ptavg = (m1.pt + m2.pt)/2;
  etaavg = (m1.eta + m2.eta)/2;
  phiavg = (m1.phi + m2.phi)/2;

  TLorentzVector M1, M2, M3;
  M1.SetPtEtaPhiM(m1.pt,m1.eta,m1.phi, 0.105658);
  M2.SetPtEtaPhiM(m2.pt,m2.eta,m2.phi, 0.105658);
  M3=M1+M2;
  minv     = M3.M();
  pair_pt  = M3.Pt();
  pair_eta = M3.Eta();
  pair_phi = M3.Phi();
  pair_y   = M3.Rapidity();

  asym = (m1.pt - m2.pt) / (m1.pt + m2.pt);
  acop = (PI - fabs(dphi)) / PI;
  if (asym < 0 || acop < 0){
    std::cout<<"Error:: Both asymmetry & acoplanarity need to be positive, quitting"<<std::endl;
    throw std::exception();
  }

  same_sign = (m1.charge == m2.charge);
  avg_centrality = (m1.ev_centrality + m2.ev_centrality)/2;
}

void MuonPair::Update(){
  Sort();
  PairValueCalc();
  // std::cout << m1.ev_num << ", parent class Update method called" << std::endl;
}

void MuonPair::UpdateShort(){
  Sort();
  double PI=acos(-1.0);

  pt_lead = m1.pt;
  pair_dPoverP = sqrt(m1.dP_overP*m1.dP_overP + m2.dP_overP*m2.dP_overP);
  dpt = m1.pt - m2.pt;
  dphi = m1.phi - m2.phi;
  dphi = atan2(sin(dphi),cos(dphi));//fold dphi to [-pi,pi]
  deta = m1.eta - m2.eta;
  dr   = sqrt(dphi*dphi + deta*deta);
  ptavg = (m1.pt + m2.pt)/2;
  etaavg = (m1.eta + m2.eta)/2;
  phiavg = (m1.phi + m2.phi)/2;

  same_sign = (m1.charge == m2.charge);
  avg_centrality = (m1.ev_centrality + m2.ev_centrality)/2;
}

#endif

