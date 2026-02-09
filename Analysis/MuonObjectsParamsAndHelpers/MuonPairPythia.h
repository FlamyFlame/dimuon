#pragma once
#include "MuonPairPbPb.h"
#include "MuonPairMC.h"

template <class Derived>
struct PairPythiaExtras {
    double crossx{};
    double effcy{};
};

template <class Derived>
struct PairPythiaTruthExtras {

    // parent groups: {direct b, c from b, c not from b, strange & light hadrons, photons}; if others: print out
    int   m1_parent_group{};
    int   m2_parent_group{};
    
    bool  from_same_b{};
    
    bool  from_same_ancestors{}; // same ancestors at the (latest non-HQ) partonic level
    
    bool  from_same_resonance{};
    bool  resonance_contaminated{};
    
    bool  Reco_resonance_or_reso_contam_tagged_old{}; // pairs tagged as either from resonances or contaminated by resonances using old Reco resonance cuts
    bool  Reco_resonance_or_reso_contam_tagged_new{}; // pairs tagged as either from resonances or contaminated by resonances using new Reco resonance cuts
    
    int   m1_hard_scatt_category{};
    int   m2_hard_scatt_category{};
    
    int   muon_pair_flavor_category{};
    int   muon_pair_origin_category{};
    
    bool  pair_origin_analysis_skipped{};
    
    bool  m1_from_pdf{};
    bool  m2_from_pdf{};
    
    double QHard{};
    double pTHat{};
    double mHat{};
    
    double Qsplit{}; // IMPORTANT: should only use for HF muon pairs from gluon/photon splitting
    double mHard_relevant{}; // IMPORTANT: should only use for HF muon pairs that are either from the same hard scattering or where only one has gone through a hard scattering
                           // If the relevant hard scattering is the hardest scattering, should be compared with mHat NOT QHard or pTHat
    // double mQQ; // IMPORTANT: should only use for HF muon pairs that are either from the same hard scattering or where only one has gone through a hard scattering
    // double dRQQ; // IMPORTANT: should only use for HF muon pairs that are either from the same hard scattering or where only one has gone through a hard scattering
    // TLorentzVector vg, vQ1, vQ2, vout1, vout2;
    // int wrong_4vec_mode_012;
    // // we'll use vQ1, vQ2 for all "same-partonic-origin" events
    // // vg for gluon splitting
    // // vout1, vout2 for (1) cross check on FC (2) same hard scattering
    
    std::vector<float> m1_last_b_hadron_prt_pt_eta_phi_m;
    std::vector<float> m2_last_b_hadron_prt_pt_eta_phi_m;
    std::vector<float> m1_last_hf_hadron_prt_pt_eta_phi_m;
    std::vector<float> m2_last_hf_hadron_prt_pt_eta_phi_m;
    std::vector<float> m1_first_hq_ancestor_pt_eta_phi_m;
    std::vector<float> m2_first_hq_ancestor_pt_eta_phi_m;
};

struct MuonPairPythiaTruth
  : MuonPairBaseT<MuonPairPythiaTruth, MuonPythia>
  , PairMCTruthExtras<MuonPairPythiaTruth>
  , PairPythiaExtras<MuonPairPythiaTruth>
  , PairPythiaTruthExtras<MuonPairPythiaTruth>
{};

struct MuonPairPythiaFullSimNoTruth
  : MuonPairBaseT<MuonPairPythiaFullSimNoTruth, MuonPythiaFullSimNoTruth>
  , PairPythiaExtras<MuonPairPythiaFullSimNoTruth>
  , PairRecoExtras<MuonPairPythiaFullSimNoTruth>
  , PairFullSimExtras<MuonPairPythiaFullSimNoTruth>
{};

struct MuonPairPythiaFullSimWTruth
  : MuonPairBaseT<MuonPairPythiaFullSimWTruth, MuonPythiaFullSimWTruth>
  , PairMCTruthExtras<MuonPairPythiaFullSimWTruth>
  , PairPythiaExtras<MuonPairPythiaFullSimWTruth>
  , PairPythiaTruthExtras<MuonPairPythiaFullSimWTruth>
  , PairRecoExtras<MuonPairPythiaFullSimWTruth>
  , PairFullSimExtras<MuonPairPythiaFullSimWTruth>
{};

struct MuonPairPythiaFullSimOverlayNoTruth
  : MuonPairBaseT<MuonPairPythiaFullSimOverlayNoTruth, MuonPythiaFullSimOverlayNoTruth>
  , PairPythiaExtras<MuonPairPythiaFullSimOverlayNoTruth>
  , PairRecoExtras<MuonPairPythiaFullSimOverlayNoTruth>
  , PairFullSimExtras<MuonPairPythiaFullSimOverlayNoTruth>
  , PairPbPbExtras<MuonPairPythiaFullSimOverlayNoTruth>
{
    void PairValueCalcHook() {
        this->PairPbPbExtras<MuonPairPythiaFullSimOverlayNoTruth>::PairValueCalcPbPb(); // compute avg_centrality
    }
};

struct MuonPairPythiaFullSimOverlayWTruth
  : MuonPairBaseT<MuonPairPythiaFullSimOverlayWTruth, MuonPythiaFullSimOverlayWTruth>
  , PairMCTruthExtras<MuonPairPythiaFullSimOverlayWTruth>
  , PairPythiaExtras<MuonPairPythiaFullSimOverlayWTruth>
  , PairPythiaTruthExtras<MuonPairPythiaFullSimOverlayWTruth>
  , PairRecoExtras<MuonPairPythiaFullSimOverlayWTruth>
  , PairFullSimExtras<MuonPairPythiaFullSimOverlayWTruth>
  , PairPbPbExtras<MuonPairPythiaFullSimOverlayWTruth>
{
    void PairValueCalcHook() {
        this->PairPbPbExtras<MuonPairPythiaFullSimOverlayWTruth>::PairValueCalcPbPb(); // compute avg_centrality
    }
};
