#ifndef MuonPairMC_h
#define MuonPairMC_h

#include "/usatlas/u/yuhanguo/workarea/dimuon_codes/MuonPair.h"
#include "vector"
#include <iostream>

class MuonPairMC : public MuonPair{ 
// define child class for the case of MC muon pairs
// with the additional attributes of parent information & helper functions to work with them
public:
  // float weight;
  float crossx;
  
  // bool m1_c_tag;
  // bool m2_c_tag;
  // bool m1_osc;
  // bool m2_osc;
  // int m1_earliest_parent_id;
  // int m2_earliest_parent_id;

  // float m1_closest_hadron_prt_pt;
  // float m2_closest_hadron_prt_pt;
  // float m1_furthest_hadron_prt_pt;
  // float m2_furthest_hadron_prt_pt;
  // float m1_hq_ancestor_pt;
  // float m2_hq_ancestor_pt;

  // parent groups: {direct b, c from b, c not from b, strange & light hadrons, photons};
  // if others: print out
  int m1_parent_group;
  int m2_parent_group;

  bool from_same_b; // only useful for the b-bbar sample, opposite sign

  MuonPairMC();
  ~MuonPairMC(){}
};

MuonPairMC::MuonPairMC(){
	// std::cout << "Constructor for MuonPairMC is called.";
}




#endif

