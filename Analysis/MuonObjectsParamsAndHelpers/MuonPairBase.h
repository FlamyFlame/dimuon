#pragma once
#include "Muon.h"
#include "TLorentzVector.h"
#include "RtypesCore.h"
#include <iostream>
#include <cstring>
#include <cmath>
#include <cassert>
#include <algorithm>   // for std::swap
#include <stdexcept>

template <class Derived, class MuonT>
struct MuonPairBaseT {
  	using muon_t = MuonT;
  	MuonT m1, m2;

  	// base/common pair quantities
	double weight{1.};
  	float pt_lead{};
	float pair_pt{}, pair_eta{}, pair_phi{}, pair_y{};
	float dpt{}, deta{}, dphi{}, dr{}; 
	float ptavg{}, etaavg{}, phiavg{};
  	float minv{};
  	float asym{}; // asymmetry := (pT_lead - pT_sublead) / (pT_lead + pT_sublead)
  	float acop{}; // acoplanarity := (pi - |Dphi|) / pi
  	bool  same_sign{};

  	Derived& self() { return static_cast<Derived&>(*this); }
  	const Derived& self() const { return static_cast<const Derived&>(*this); }

  	void Sort() {
        if (m2.pt > m1.pt) {
            std::swap(m1, m2);
        }
    }

  	void PairValueCalc() {
      	PairValueCalcBase();      // always do base part
    	self().PairValueCalcHook(); // derived hook (provided or defaulted)
  	}

    void Update(){
        Sort();
        PairValueCalc();
    }

protected:
  	void PairValueCalcBase();
};

template <class Derived>
struct PairCalcHookDefault {
  	void PairValueCalcHook() {}  // no-op by default
};


template <class Derived, class MuonT>
void MuonPairBaseT<Derived, MuonT>::PairValueCalcBase(){
    double PI=acos(-1.0);

    pt_lead = m1.pt;
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

    same_sign = (m1.charge == m2.charge);
}
