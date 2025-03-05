#ifndef Muon_h
#define Muon_h
#include <vector>

class Muon{
public:
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

   	Muon(){}
   	~Muon(){}
};

#endif