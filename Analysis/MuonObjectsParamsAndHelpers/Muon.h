#pragma once
#include <vector>

// -------------------- Muon helper structs (base & extras) --------------------
struct MuonBase {
  	float pt{-1e6}, eta{1e6}, phi{1e6};
  	int   ind{-1}, ev_num{-1}, charge{0};
};

struct MuonData : MuonBase{
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

struct MuonPbPb: MuonData, MuonPbPbExtra{};

struct MuonMCExtra {
    float truth_pt{-1e6}, truth_eta{1e6}, truth_phi{1e6};
    int   truth_bar{-1}, truth_ch{0};
};

struct MuonPythiaExtra {
    // pythia-only fields (placeholder)
};

struct MuonPowhegExtra {
    // powheg-only fields (placeholder)
};

struct MuonPythiaTruth : MuonBase, MuonMCExtra, MuonPythiaExtra {};
struct MuonPowhegTruth : MuonBase, MuonMCExtra, MuonPowhegExtra {};

using MuonPythia = MuonPythiaTruth;
using MuonPowheg = MuonPowhegTruth;

struct MuonPythiaFullSim : MuonData, MuonMCExtra, MuonPythiaExtra {};
struct MuonPowhegFullSim : MuonData, MuonMCExtra, MuonPowhegExtra {};

struct MuonPythiaFullSimOverlay : MuonData, MuonMCExtra, MuonPythiaExtra, MuonPbPbExtra {};
struct MuonPowhegFullSimOverlay : MuonData, MuonMCExtra, MuonPowhegExtra, MuonPbPbExtra {};

