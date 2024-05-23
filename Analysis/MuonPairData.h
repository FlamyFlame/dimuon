#ifndef MuonPairData_h
#define MuonPairData_h

#include "/usatlas/u/yuhanguo/workarea/dimuon_codes/MuonPair.h"
#include "vector"
#include <iostream>

class MuonPairData : public MuonPair{ 
// define child class for the case of MC muon pairs
// with the additional attributes of parent information & helper functions to work with them
public:
  float m1_trk_pt;
  float m1_trk_eta;
  float m1_trk_phi;
  float m2_trk_pt;
  float m2_trk_eta;
  float m2_trk_phi;

  MuonPairData(){}
  ~MuonPairData(){}
  void Sort();
  void Update();
};

void MuonPairData::Sort(){
  ///sort pt, eta, phi by pt

  float temppt, tempeta, tempphi, tempd0, tempz0, tempdP_overP;
  int tempind, tempcharge, tempquality;
  int temp_trkpt, temp_trketa, temp_trkphi;

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
    temp_trkpt = m1_trk_pt;
    temp_trketa = m1_trk_eta;
    temp_trkphi = m1_trk_phi;
  
    m1.pt = m2.pt;
    m1.eta = m2.eta;
    m1.phi = m2.phi;
    m1.d0 = m2.d0;
    m1.z0 = m2.z0;
    m1.ind = m2.ind;
    m1.charge = m2.charge;
    m1.quality = m2.quality;
    m1.dP_overP = m2.dP_overP;
    m1_trk_pt = m2_trk_pt;
    m1_trk_eta = m2_trk_eta;
    m1_trk_phi = m2_trk_phi;
  
    m2.pt = temppt;
    m2.eta = tempeta;
    m2.phi = tempphi;
    m2.d0 = tempd0;
    m2.z0 = tempz0;
    m2.ind = tempind;
    m2.charge = tempcharge;
    m2.quality = tempquality;
    m2.dP_overP = tempdP_overP;
    m1_trk_pt = temp_trkpt;
    m1_trk_eta = temp_trketa;
    m1_trk_phi = temp_trkphi;
  }
  // std::cout << m1.ev_num << ", child class sorting method called" << std::endl;
}

void MuonPairData::Update(){
  Sort(); // use the overriden version that also sorts the track pt, eta, phi
  PairValueCalc(); // use the parent-class implementation
  // std::cout << m1.ev_num << ", child class update method called" << std::endl;
}

#endif

