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
    bool    pass_tight{};

    float   trk_pt{-1e6}, trk_eta{1e6}, trk_phi{1e6};
    int     trk_charge{};
};

struct MuonPbPbExtra{
   	int ev_centrality{-1};
   	float ev_FCal_Et{-1e6};
};

struct MuonPP: MuonBase, MuonRecoExtra{};
struct MuonPbPb: MuonBase, MuonRecoExtra, MuonPbPbExtra{};

struct MuonMCTruthKinExtra {
    float truth_pt{-1e6}, truth_eta{1e6}, truth_phi{1e6};
    int   truth_bar{-1}, truth_charge{0};
};

struct MuonFullsimExtra {
    bool pass_d0_z0{}, pass_medium{};
    bool reco_match{}; // matched with a reco muon with prob > 0.5
};

struct MuonPythiaExtra {
    // pythia-only fields (placeholder)
};

struct MuonPowhegExtra {
    // powheg-only fields (placeholder)
};

struct MuonPythiaTruth : MuonBase, MuonMCTruthKinExtra, MuonPythiaExtra {};
struct MuonPowhegTruth : MuonBase, MuonMCTruthKinExtra, MuonPowhegExtra {};

struct MuonPythiaFullSimNoTruth : MuonBase, MuonRecoExtra, MuonMCTruthKinExtra, MuonPythiaExtra {};
struct MuonPowhegFullSimNoTruth : MuonBase, MuonRecoExtra, MuonMCTruthKinExtra, MuonPowhegExtra {};

struct MuonPythiaFullSimWTruth : MuonBase, MuonRecoExtra, MuonMCTruthKinExtra, MuonPythiaExtra {};
struct MuonPowhegFullSimWTruth : MuonBase, MuonRecoExtra, MuonMCTruthKinExtra, MuonPowhegExtra {};

struct MuonPythiaFullSimOverlayNoTruth : MuonBase, MuonRecoExtra, MuonMCTruthKinExtra, MuonPythiaExtra, MuonPbPbExtra {};
struct MuonPowhegFullSimOverlayNoTruth : MuonBase, MuonRecoExtra, MuonMCTruthKinExtra, MuonPowhegExtra, MuonPbPbExtra {};

struct MuonPythiaFullSimOverlayWTruth : MuonBase, MuonRecoExtra, MuonMCTruthKinExtra, MuonPythiaExtra, MuonPbPbExtra {};
struct MuonPowhegFullSimOverlayWTruth : MuonBase, MuonRecoExtra, MuonMCTruthKinExtra, MuonPowhegExtra, MuonPbPbExtra {};

