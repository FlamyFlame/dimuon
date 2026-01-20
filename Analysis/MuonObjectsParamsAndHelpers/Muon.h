#pragma once
#include <vector>

// -------------------- Muon base --------------------
struct MuonBase {
  	int   ind{-1}, ev_num{-1};
};

// -------------------- Muon reco --------------------
struct MuonRecoExtra{
    float pt{-1e6}, eta{1e6}, phi{1e6}; // reco kinematics
    int   charge{0};

   	float   dP_overP{1e6}, z0{1e6}, d0{1e6};
   	int     quality{-1};
	bool    passmu4{};
    bool    passTight{};

    float   trk_pt{-1e6}, trk_eta{1e6}, trk_phi{1e6};
    int     trk_charge{};
};

struct MuonPbPbExtra{
   	int ev_centrality{-1};
   	float ev_FCal_Et{-1e6};
};

struct MuonPP: MuonBase, MuonRecoExtra{};
struct MuonPbPb: MuonBase, MuonRecoExtra, MuonPbPbExtra{};

struct MuonMCTruthExtra {
    float truth_pt{-1e6}, truth_eta{1e6}, truth_phi{1e6};
    int   truth_bar{-1}, truth_charge{0};
};

struct MuonPythiaExtra {
    // pythia-only fields (placeholder)
};

struct MuonPowhegExtra {
    // powheg-only fields (placeholder)
};

struct MuonPythiaTruth : MuonBase, MuonMCTruthExtra, MuonPythiaExtra {};
struct MuonPowhegTruth : MuonBase, MuonMCTruthExtra, MuonPowhegExtra {};

struct MuonPythiaFullSimNoTruth : MuonBase, MuonRecoExtra, MuonPythiaExtra {};
struct MuonPowhegFullSimNoTruth : MuonBase, MuonRecoExtra, MuonPowhegExtra {};

struct MuonPythiaFullSimWTruth : MuonBase, MuonRecoExtra, MuonMCTruthExtra, MuonPythiaExtra {};
struct MuonPowhegFullSimWTruth : MuonBase, MuonRecoExtra, MuonMCTruthExtra, MuonPowhegExtra {};

struct MuonPythiaFullSimOverlayNoTruth : MuonBase, MuonRecoExtra, MuonPythiaExtra, MuonPbPbExtra {};
struct MuonPowhegFullSimOverlayNoTruth : MuonBase, MuonRecoExtra, MuonPowhegExtra, MuonPbPbExtra {};

struct MuonPythiaFullSimOverlayWTruth : MuonBase, MuonRecoExtra, MuonMCTruthExtra, MuonPythiaExtra, MuonPbPbExtra {};
struct MuonPowhegFullSimOverlayWTruth : MuonBase, MuonRecoExtra, MuonMCTruthExtra, MuonPowhegExtra, MuonPbPbExtra {};

