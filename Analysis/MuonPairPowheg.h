#ifndef MuonPairPowheg_h
#define MuonPairPowheg_h

#include "MuonPair.h"
#include "vector"
#include <iostream>

class MuonPairPowheg : public MuonPair{ 
// define child class for the case of MC muon pairs
// with the additional attributes of parent information & helper functions to work with them
public:
  // float weight;
  float crossx;
  

  // parent groups: {direct b, c from b, c not from b, strange & light hadrons, photons};
  // if others: print out
  int m1_parent_group;
  int m2_parent_group;

  bool from_same_b; // only useful for the b-bbar sample, opposite sign

  bool  from_same_ancestors;
  bool  both_from_b;
  // bool  one_from_b_one_from_c;
  bool  both_from_c;

  int   m1_ancestor_category;
  int   m2_ancestor_category;

  float Q;
  float mQQ; // mHat for the relevant hard scattering
  float mHard_relevant; // mHat for the relevant hard scattering
  // float s_cm; // mHat for the relevant hard scattering

  std::vector<float> m1_last_b_hadron_prt_pt_eta_phi_m;
  std::vector<float> m2_last_b_hadron_prt_pt_eta_phi_m;
  std::vector<float> m1_last_hf_hadron_prt_pt_eta_phi_m;
  std::vector<float> m2_last_hf_hadron_prt_pt_eta_phi_m;
  std::vector<float> m1_first_hq_ancestor_pt_eta_phi_m;
  std::vector<float> m2_first_hq_ancestor_pt_eta_phi_m;


  MuonPairPowheg();
  ~MuonPairPowheg(){}
};

MuonPairPowheg::MuonPairPowheg(){
	// std::cout << "Constructor for MuonPairPowheg is called.";
}


#endif