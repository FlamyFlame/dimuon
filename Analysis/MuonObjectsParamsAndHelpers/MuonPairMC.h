#pragma once
#include "MuonPairBase.h"

template <class Derived>
struct PairMCTruthExtras {
    float truth_pt_lead{};
    float truth_pair_pt{}, truth_pair_eta{}, truth_pair_phi{}, truth_pair_y{};
    float truth_dpt{}, truth_deta{}, truth_dphi{}, truth_dr{}; 
    float truth_ptavg{}, truth_etaavg{}, truth_phiavg{};
    bool  truth_same_sign{};

    float truth_minv{}, truth_asym{}, truth_acop{};

    void SortTruth() {
        auto& d = static_cast<Derived&>(*this);
        if (d.m2.truth_pt > d.m1.truth_pt) std::swap(d.m1, d.m2);
    }

    void PairValueCalcTruth(){
        auto& d = static_cast<Derived&>(*this);

        double PI=acos(-1.0);

        truth_pt_lead = d.m1.truth_pt;
        truth_dpt = d.m1.truth_pt - d.m2.truth_pt;
        truth_dphi = d.m1.truth_phi - d.m2.truth_phi;
        truth_dphi = atan2(sin(truth_dphi),cos(truth_dphi));//fold truth_dphi to [-pi,pi]
        truth_deta = d.m1.truth_eta - d.m2.truth_eta;
        truth_dr   = sqrt(truth_dphi * truth_dphi + truth_deta * truth_deta);
        truth_ptavg = (d.m1.truth_pt + d.m2.truth_pt)/2;
        truth_etaavg = (d.m1.truth_eta + d.m2.truth_eta)/2;
        truth_phiavg = (d.m1.truth_phi + d.m2.truth_phi)/2;

        TLorentzVector M1, M2, M3;
        M1.SetPtEtaPhiM(d.m1.truth_pt,d.m1.truth_eta,d.m1.truth_phi, 0.105658);
        M2.SetPtEtaPhiM(d.m2.truth_pt,d.m2.truth_eta,d.m2.truth_phi, 0.105658);
        M3=M1+M2;
        truth_minv     = M3.M();
        truth_pair_pt  = M3.Pt();
        truth_pair_eta = M3.Eta();
        truth_pair_phi = M3.Phi();
        truth_pair_y   = M3.Rapidity();

        truth_asym = (d.m1.truth_pt - d.m2.truth_pt) / (d.m1.truth_pt + d.m2.truth_pt);
        truth_acop = (PI - fabs(truth_dphi)) / PI;

        truth_same_sign = (d.m1.truth_charge == d.m2.truth_charge);
    }
};


