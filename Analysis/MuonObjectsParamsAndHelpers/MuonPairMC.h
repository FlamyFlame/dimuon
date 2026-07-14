#pragma once
#include "MuonPairBase.h"

template <class Derived>
struct PairMCTruthKinExtras {
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

template <class Derived>
struct PairFullSimExtras {
    bool pair_pass_medium{}; // default medium pair selection WITHOUT resonance cuts
    bool pair_pass_medium_and_resonance{};
    bool pair_pass_tight_and_resonance{};
    bool pair_pass_resonance_reco{};
    bool pair_pass_resonance_truth{};
};

// Pair-level trigger decisions propagated from the raw-NTUP trigger blocks
// (filled only in store_mc_trigger mode; default false otherwise). The per-leg mu4
// decisions live on the muons (m1.passmu4 / m2.passmu4, MuonRecoExtra).
// NOTE: keep >= 2 data members — ROOT does not member-split a single-member base
// class into named leaves (the branch gets one leaf named after the base, which
// RDataFrame cannot address by member name).
template <class Derived>
struct PairMCTrigExtras {
    bool pass2mu4{};    // order-insensitive dimuon_b_HLT_2mu4_L12MU3V match (mindR 0.02, data-mirror)
    bool ev_pass_mu4{}; // event-level b_HLT_mu4_L1MU3V decision (diagnostics: mu4 superset checks)
    bool ev_pass_2mu4{};// event-level b_HLT_2mu4_L12MU3V decision (diagnostics: 2mu4 vs per-pair match)
};
