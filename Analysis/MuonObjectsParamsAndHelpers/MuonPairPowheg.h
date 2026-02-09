#pragma once
#include "MuonPairPbPb.h"
#include "MuonPairMC.h"

template <class Derived>
struct PairPowhegExtras {
    float crossx;
};

template <class Derived>
struct PairPowhegTruthExtras {
    float Q{-1000.};

    // parent groups: {direct b, c from b, c not from b, strange & light hadrons, photons}; if others: print out
    int m1_parent_group{-1};
    int m2_parent_group{-1};

    bool from_same_b{false}; // only useful for the b-bbar sample, opposite sign

    bool  from_same_ancestors{false};
    bool  both_from_b{false};
    // bool  one_from_b_one_from_c;
    bool  both_from_c{false};

    int   m1_ancestor_category{-10};
    int   m2_ancestor_category{-10};

    float mQQ{-10.}; // mHat for the relevant hard scattering
    float mHard_relevant{-10.}; // mHat for the relevant hard scattering
    // float s_cm; // mHat for the relevant hard scattering

    std::vector<float> m1_last_b_hadron_prt_pt_eta_phi_m;
    std::vector<float> m2_last_b_hadron_prt_pt_eta_phi_m;
    std::vector<float> m1_last_hf_hadron_prt_pt_eta_phi_m;
    std::vector<float> m2_last_hf_hadron_prt_pt_eta_phi_m;
    std::vector<float> m1_first_hq_ancestor_pt_eta_phi_m;
    std::vector<float> m2_first_hq_ancestor_pt_eta_phi_m;
};

struct MuonPairPowhegTruth
  : MuonPairBaseT<MuonPairPowhegTruth, MuonPowhegTruth>
  , PairMCTruthKinExtras<MuonPairPowhegTruth>
  , PairPowhegExtras<MuonPairPowhegTruth>
  , PairPowhegTruthExtras<MuonPairPowhegTruth>
{};

struct MuonPairPowhegFullSimNoTruth
  : MuonPairBaseT<MuonPairPowhegFullSimNoTruth, MuonPowhegFullSimNoTruth>
  , PairMCTruthKinExtras<MuonPairPowhegFullSimNoTruth>
  , PairPowhegExtras<MuonPairPowhegFullSimNoTruth>
  , PairRecoExtras<MuonPairPowhegFullSimNoTruth>
  , PairFullSimExtras<MuonPairPowhegFullSimNoTruth>
{};

struct MuonPairPowhegFullSimWTruth
  : MuonPairBaseT<MuonPairPowhegFullSimWTruth, MuonPowhegFullSimWTruth>
  , PairMCTruthKinExtras<MuonPairPowhegFullSimWTruth>
  , PairPowhegExtras<MuonPairPowhegFullSimWTruth>
  , PairPowhegTruthExtras<MuonPairPowhegFullSimWTruth>
  , PairRecoExtras<MuonPairPowhegFullSimWTruth>
  , PairFullSimExtras<MuonPairPowhegFullSimWTruth>
{};

struct MuonPairPowhegFullSimOverlayNoTruth
  : MuonPairBaseT<MuonPairPowhegFullSimOverlayNoTruth, MuonPowhegFullSimOverlayNoTruth>
  , PairMCTruthKinExtras<MuonPairPowhegFullSimOverlayNoTruth>
  , PairPowhegExtras<MuonPairPowhegFullSimOverlayNoTruth>
  , PairRecoExtras<MuonPairPowhegFullSimOverlayNoTruth>
  , PairFullSimExtras<MuonPairPowhegFullSimOverlayNoTruth>
  , PairPbPbExtras<MuonPairPowhegFullSimOverlayNoTruth>
{
    void PairValueCalcHook() {
        this->PairPbPbExtras<MuonPairPowhegFullSimOverlayNoTruth>::PairValueCalcPbPb(); // compute avg_centrality
    }
};

struct MuonPairPowhegFullSimOverlayWTruth
  : MuonPairBaseT<MuonPairPowhegFullSimOverlayWTruth, MuonPowhegFullSimOverlayWTruth>
  , PairMCTruthKinExtras<MuonPairPowhegFullSimOverlayWTruth>
  , PairPowhegExtras<MuonPairPowhegFullSimOverlayWTruth>
  , PairPowhegTruthExtras<MuonPairPowhegFullSimOverlayWTruth>
  , PairRecoExtras<MuonPairPowhegFullSimOverlayWTruth>
  , PairFullSimExtras<MuonPairPowhegFullSimOverlayWTruth>
  , PairPbPbExtras<MuonPairPowhegFullSimOverlayWTruth>
{
    void PairValueCalcHook() {
        this->PairPbPbExtras<MuonPairPowhegFullSimOverlayWTruth>::PairValueCalcPbPb();
    }
};
