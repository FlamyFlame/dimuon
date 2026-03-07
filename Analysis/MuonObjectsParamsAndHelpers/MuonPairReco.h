#pragma once
#include "MuonPairBase.h"


template <class Derived>
struct PairDataExtras {
    UInt_t  run_number{}, lb{}, bcid{};
    
    bool    passmu4mu4noL1;
    bool    passmu4noL1; // pair-level: true if either muon in the pair passes the mu4noL1 (unseeded) leg
    bool    pass2mu4;
    bool    passSeparated; // separated enough for 2mu4 & mu4_mu4noL1 trigger efficiencies to be factorizable into single-muon parts; for now dR > 0.8
    bool    passSeparatedDeta; // separated enough for 2mu4 & mu4_mu4noL1 trigger efficiencies to be factorizable into single-muon parts; for now deta > 0.8
};

template <class Derived>
struct PairRecoExtras {
  	float pt_lead{};
    float pair_pt{}, pair_eta{}, pair_phi{}, pair_y{};
    float dpt{}, deta{}, dphi{}, dr{}; 
    float ptavg{}, etaavg{}, phiavg{};
    float minv{};
    float asym{}; // asymmetry := (pT_lead - pT_sublead) / (pT_lead + pT_sublead)
    float acop{}; // acoplanarity := (pi - |Dphi|) / pi
    bool  same_sign{};

    float pair_dPoverP{};

    bool    pair_pass_tight;

    void SortReco() {
        auto& d = static_cast<Derived&>(*this);
        if (d.m2.pt > d.m1.pt) std::swap(d.m1, d.m2);
    }

	void PairValueCalcReco() {
      	auto& d = static_cast<Derived&>(*this);

        double PI=acos(-1.0);

        pt_lead = d.m1.pt;
        dpt = d.m1.pt - d.m2.pt;
        dphi = d.m1.phi-d.m2.phi;
        dphi = atan2(sin(dphi),cos(dphi));//fold dphi to [-pi,pi]
        deta = d.m1.eta-d.m2.eta;
        dr   = sqrt(dphi*dphi + deta*deta);
        ptavg = (d.m1.pt + d.m2.pt)/2;
        etaavg = (d.m1.eta + d.m2.eta)/2;
        phiavg = (d.m1.phi + d.m2.phi)/2;

        TLorentzVector M1, M2, M3;
        M1.SetPtEtaPhiM(d.m1.pt,d.m1.eta,d.m1.phi, 0.105658);
        M2.SetPtEtaPhiM(d.m2.pt,d.m2.eta,d.m2.phi, 0.105658);
        M3=M1+M2;
        minv     = M3.M();
        pair_pt  = M3.Pt();
        pair_eta = M3.Eta();
        pair_phi = M3.Phi();
        pair_y   = M3.Rapidity();

        asym = (d.m1.pt - d.m2.pt) / (d.m1.pt + d.m2.pt);
        acop = (PI - fabs(dphi)) / PI;

        same_sign = (d.m1.charge == d.m2.charge);
    	pair_dPoverP = std::sqrt(d.m1.dP_overP*d.m1.dP_overP + d.m2.dP_overP*d.m2.dP_overP);
  	}
};

struct MuonPairPP
  : MuonPairBaseT<MuonPairPP, MuonPP>
  , PairRecoExtras<MuonPairPP>
  , PairDataExtras<MuonPairPP>
{};
