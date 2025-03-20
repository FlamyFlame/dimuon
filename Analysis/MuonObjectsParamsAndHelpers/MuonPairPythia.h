#ifndef MuonPairPythia_h
#define MuonPairPythia_h

#include "MuonPair.h"
#include "vector"
#include <iostream>

class MuonPairPythia : public MuonPair{ 
// define child class for the case of MC muon pairs
// with the additional attributes of parent information & helper functions to work with them
public:
  // float weight;
  double crossx;
  double effcy;
  
  // parent groups: {direct b, c from b, c not from b, strange & light hadrons, photons};
  // if others: print out
  int   m1_parent_group;
  int   m2_parent_group;

  bool  from_same_b;
  bool  from_drell_yan;
  bool  from_photo2dimuon;

  // int hard_scatt_type;

  bool  from_same_ancestors;
  bool  both_from_b;
  bool  one_from_b_one_from_c;
  bool  both_from_c;
  bool  from_same_resonance;
  bool  resonance_contaminated;

  int   m1_hard_scatt_category;
  int   m2_hard_scatt_category;
  int   muon_pair_origin_category;

  bool m1_from_pdf;
  bool m2_from_pdf;

  double QHard;
  double pTHat;
  double mHat;

  double Qsplit; // IMPORTANT: should only use for HF muon pairs from gluon/photon splitting
  double mHard_relevant; // IMPORTANT: should only use for HF muon pairs that are either from the same hard scattering or where only one has gone through a hard scattering
                         // If the relevant hard scattering is the hardest scattering, should be compared with mHat NOT QHard or pTHat
  // double mQQ; // IMPORTANT: should only use for HF muon pairs that are either from the same hard scattering or where only one has gone through a hard scattering
  // double dRQQ; // IMPORTANT: should only use for HF muon pairs that are either from the same hard scattering or where only one has gone through a hard scattering
  // TLorentzVector vg, vQ1, vQ2, vout1, vout2;
  // int wrong_4vec_mode_012;
  // // we'll use vQ1, vQ2 for all "same-partonic-origin" events
  // // vg for gluon splitting
  // // vout1, vout2 for (1) cross check on FC (2) same hard scattering
  
  std::vector<float> m1_last_b_hadron_prt_pt_eta_phi_m;
  std::vector<float> m2_last_b_hadron_prt_pt_eta_phi_m;
  std::vector<float> m1_last_hf_hadron_prt_pt_eta_phi_m;
  std::vector<float> m2_last_hf_hadron_prt_pt_eta_phi_m;
  std::vector<float> m1_first_hq_ancestor_pt_eta_phi_m;
  std::vector<float> m2_first_hq_ancestor_pt_eta_phi_m;

  // bool from_same_b; // only useful for the b-bbar sample, opposite sign

  MuonPairPythia();
  ~MuonPairPythia(){}
};

MuonPairPythia::MuonPairPythia(){
	// std::cout << "Constructor for MuonPairPythia is called.";
}


#endif
