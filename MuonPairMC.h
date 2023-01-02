#ifndef MuonPairMC_h
#define MuonPairMC_h

#include "MuonPair.h"
#include "vector"
#include <iostream>

class MuonPairMC : public MuonPair{ 
// define child class for the case of MC muon pairs
// with the additional attributes of parent information & helper functions to work with them
public:
  // std::vector<int> m1_parent_ids;
  // std::vector<int> m1_parent_barcodes;
  // std::vector<int> m2_parent_ids;
  // std::vector<int> m2_parent_barcodes;
  int m1_parent_id;
  int m1_parent_barcode;
  int m2_parent_id;
  int m2_parent_barcode;

  bool m1_c_tag;
  bool m2_c_tag;
  // std::vector<int> m1_first_non_c_parent_ids;
  // std::vector<int> m1_first_non_c_parent_barcodes;
  // std::vector<int> m2_first_non_c_parent_ids;
  // std::vector<int> m2_first_non_c_parent_barcodes;
  int m1_earliest_parent_id;
  // int m1_first_non_c_parent_barcode;
  int m2_earliest_parent_id;
  // int m2_first_non_c_parent_barcode;

  // parent groups: {direct b, c from b, c not from b, strange & light hadrons, photons};
  // if others: print out
  int m1_parent_group;
  int m2_parent_group;


  MuonPairMC();
  ~MuonPairMC(){}
};

MuonPairMC::MuonPairMC(){
	std::cout << "Constructor for MuonPairMC is called.";
}




#endif

