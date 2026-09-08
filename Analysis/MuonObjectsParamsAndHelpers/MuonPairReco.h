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

// pp-ONLY pair extras. Deliberately NOT in PairDataExtras, which Pb+Pb also inherits: Pb+Pb
// rejects pile-up at event level and its pairs are primary-vertex pairs by construction, so it
// must neither carry nor be asked about this field (and its output tree stays unchanged).
template <class Derived>
struct PairPPExtras {
    // Index, in the skim's PrimaryVertices dump (branches vtx_x/vtx_y/vtx_z/vtx_ntrk, stored in
    // container order), of the vertex with respect to which BOTH muons of this pair pass the
    // |d0| and |z0 sin(theta)| cuts -- the best-matching such vertex (see
    // PPExtras::PassD0Z0Extra). 0 IS the primary vertex; anything > 0 is a secondary (pile-up)
    // vertex. -1 means "not set", which for a pair that reached the output tree is a bug.
    //
    // WHY pp keeps pairs off the primary vertex at all: pp24 has ~4 inelastic collisions per
    // bunch crossing and the pp luminosity is defined over ALL of them, so counting only
    // primary-vertex pairs puts an under-counted numerator over an all-collision L_int.
    // See docs/tracking/pp24_all_vertex_pairs.md.
    int matched_vtx_ind{-1};

    // Would this pair ALSO have passed the old primary-vertex-only cut? Kept so the two
    // distinct statistics can both be computed from the ntuple OUTPUT and are never confused:
    //   pairs "from a secondary vertex"        <=>  matched_vtx_ind != 0
    //   pairs ADDED by the all-vertex change   <=>  !pass_primary_vtx
    // They differ: a pair can pass w.r.t. the primary and still match a secondary vertex more
    // closely, in which case it is the former but not the latter.
    bool pass_primary_vtx{false};
};

struct MuonPairPP
  : MuonPairBaseT<MuonPairPP, MuonPP>
  , PairRecoExtras<MuonPairPP>
  , PairDataExtras<MuonPairPP>
  , PairPPExtras<MuonPairPP>
{};
