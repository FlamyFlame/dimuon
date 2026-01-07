#pragma once
#include "MuonPairBase.h"


template <class Derived>
struct PairDataExtras {
  	float pair_dPoverP{};

    UInt_t  run_number{}, lb{}, bcid{};
    
    bool    passTight;
    bool    passmu4mu4noL1;
    bool    pass2mu4;
    bool    passSeparated; // separated enough for 2mu4 & mu4_mu4noL1 trigger efficiencies to be factorizable into single-muon parts; for now dR > 0.8
    bool    passSeparatedDeta; // separated enough for 2mu4 & mu4_mu4noL1 trigger efficiencies to be factorizable into single-muon parts; for now deta > 0.8

	void PairValueCalcData() {
      	auto& d = static_cast<Derived&>(*this);
    	d.pair_dPoverP = std::sqrt(d.m1.dP_overP*d.m1.dP_overP + d.m2.dP_overP*d.m2.dP_overP);
  	}
};

struct MuonPairData
  : MuonPairBaseT<MuonPairData, MuonData>
  , PairDataExtras<MuonPairData>
  , PairCalcHookDefault<MuonPairData>
{
    void PairValueCalcHook() {
        this->PairValueCalcData();
    }
};