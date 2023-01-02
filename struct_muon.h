#ifndef struct_muon_h
#define struct_muon_h
#include <vector>

struct Muon{
	float pt;
   	float eta;
   	float phi;
   	float dP_overP;
   	float z0;
   	float d0;
   	int ind;
   	int charge;
   	int quality;
   	int ev_num;
   	int ev_centrality;
   	float ev_FCal_Et;
};

struct TruthMuon{
	float pt;
   	float eta;
   	float phi;
   	// float dP_overP;
   	int charge;
   	int barcode; //MC
   	std::vector<int> parent_ids;
   	std::vector<int> parent_barcodes;
   	std::vector<std::vector<int>> grandparent_ids;
   	std::vector<std::vector<int>> grandparent_barcodes;
};

#endif