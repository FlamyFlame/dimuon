#pragma once
#include <vector>

// -------------------- Muon helper structs (base & extras) --------------------
struct MuonBase {
  	int   ind{-1}, ev_num{-1};
};

struct MuonReco : MuonBase{
    float pt{-1e6}, eta{1e6}, phi{1e6}; // reco kinematics
    int   charge{0};

   	float   dP_overP{1e6}, z0{1e6}, d0{1e6};
   	int     quality{-1};
	bool    passmu4{};

    float   trk_pt{-1e6}, trk_eta{1e6}, trk_phi{1e6};
    int     trk_charge{};
};

struct MuonPbPbExtra{
   	int ev_centrality{-1};
   	float ev_FCal_Et{-1e6};
};

using MuonPP = MuonReco;

struct MuonPbPb: MuonReco, MuonPbPbExtra{};

struct MuonMCTruthExtra {
    float truth_pt{-1e6}, truth_eta{1e6}, truth_phi{1e6};
    int   truth_ind{-1}, truth_charge{0};
};

struct MuonPythiaExtra {
    // pythia-only fields (placeholder)
};

struct MuonPowhegExtra {
    // powheg-only fields (placeholder)
};

struct MuonPythiaTruth : MuonBase, MuonMCTruthExtra, MuonPythiaExtra {};
struct MuonPowhegTruth : MuonBase, MuonMCTruthExtra, MuonPowhegExtra {};

struct MuonPythiaFullSimNoTruth : MuonReco, MuonPythiaExtra {};
struct MuonPowhegFullSimNoTruth : MuonReco, MuonPowhegExtra {};

struct MuonPythiaFullSimWTruth : MuonReco, MuonMCTruthExtra, MuonPythiaExtra {};
struct MuonPowhegFullSimWTruth : MuonReco, MuonMCTruthExtra, MuonPowhegExtra {};

struct MuonPythiaFullSimOverlayNoTruth : MuonReco, MuonPythiaExtra, MuonPbPbExtra {};
struct MuonPowhegFullSimOverlayNoTruth : MuonReco, MuonPowhegExtra, MuonPbPbExtra {};

struct MuonPythiaFullSimOverlayWTruth : MuonReco, MuonMCTruthExtra, MuonPythiaExtra, MuonPbPbExtra {};
struct MuonPowhegFullSimOverlayWTruth : MuonReco, MuonMCTruthExtra, MuonPowhegExtra, MuonPbPbExtra {};

