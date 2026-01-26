#pragma once
#include "MuonPairPbPb.h"
#include "MuonPairMC.h"

template <class Derived>
struct PairPowhegExtras {
    float crossx;
    float Q;
};

template <class Derived>
struct PairPowhegTruthExtras {

    // parent groups: {direct b, c from b, c not from b, strange & light hadrons, photons}; if others: print out
    int m1_parent_group;
    int m2_parent_group;

    bool from_same_b; // only useful for the b-bbar sample, opposite sign

    bool  from_same_ancestors;
    bool  both_from_b;
    // bool  one_from_b_one_from_c;
    bool  both_from_c;

    int   m1_ancestor_category;
    int   m2_ancestor_category;

    float mQQ; // mHat for the relevant hard scattering
    float mHard_relevant; // mHat for the relevant hard scattering
    // float s_cm; // mHat for the relevant hard scattering

    std::vector<float> m1_last_b_hadron_prt_pt_eta_phi_m;
    std::vector<float> m2_last_b_hadron_prt_pt_eta_phi_m;
    std::vector<float> m1_last_hf_hadron_prt_pt_eta_phi_m;
    std::vector<float> m2_last_hf_hadron_prt_pt_eta_phi_m;
    std::vector<float> m1_first_hq_ancestor_pt_eta_phi_m;
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
