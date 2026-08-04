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
    bool    passmu4noL1{}; // passes the mu4noL1 (unseeded) leg of the mu4_mu4noL1 dimuon trigger
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
    float ev_weight{0.0f};
};

struct MuonFullsimExtra {
    bool pass_medium{};
    // Per-muon L1 leg for the MC trigger-efficiency L1/HLT split (mc_trigger_efficiency.md
    // round-5 change #3): true iff the offline muon matched an L1_MU3V muon RoI
    // (skim branch muon_match_L1MU3V). Filled only in store_mc_trigger mode, and only when
    // the skim carries the branch (re-skimmed NTUPs); default false otherwise. Prescale-free.
    // eff(L1) = P[pass_l1 | offline]; eff(HLT|L1) = P[passmu4 (full chain) | pass_l1].
    bool pass_l1{};
    bool reco_match{}; // matched with a reco muon with prob > 0.5
    int  reco_ind{-1}; // index of the matched reco muon in the raw-NTUP muon block (-1 = unmatched);
                       // needed to look up per-muon and per-pair trigger-matching branches
    // reco/ID efficiency scale factors (skim MuonEfficiencyScaleFactors tools, WP-dependent;
    // the skim fills them only for muons passing that WP, <=0 otherwise). Unfilled -> 1
    // (counted in the NTP printout). Filled only in store_mc_trigger mode; default 1.
    float eff_sf_medium{1.f};
    float eff_sf_tight{1.f};
    // Number of RECONSTRUCTED, TRACK-BEARING primary vertices in the event (round 7,
    // mc_trigger_efficiency.md §3.5 sanity check): count of `vtx_ntrk >= 2` entries in the
    // skim's PrimaryVertices dump. The skim stores the container unfiltered, so it always
    // ALSO contains exactly one dummy beamspot vertex with ntrk = 0, which is NOT counted.
    // n_vtx == 1 selects pile-up-free events (pp24 fullsim: ~14.5% of events; the HIJING
    // overlay reconstructs NO track-based vertex at all, so n_vtx == 0 there and the
    // requirement is inapplicable). Filled only in store_mc_trigger mode; -1 otherwise.
    int n_vtx{-1};
};

struct MuonPythiaExtra {
    // pythia-only fields (placeholder)
};

struct MuonPowhegExtra {
    // powheg-only fields (placeholder)
};

struct MuonPythiaTruth : MuonBase, MuonMCTruthKinExtra, MuonPythiaExtra {};
struct MuonPowhegTruth : MuonBase, MuonMCTruthKinExtra, MuonPowhegExtra {};

struct MuonPythiaFullSimNoTruth : MuonBase, MuonRecoExtra, MuonFullsimExtra, MuonMCTruthKinExtra, MuonPythiaExtra {};
struct MuonPowhegFullSimNoTruth : MuonBase, MuonRecoExtra, MuonFullsimExtra, MuonMCTruthKinExtra, MuonPowhegExtra {};

struct MuonPythiaFullSimWTruth : MuonBase, MuonRecoExtra, MuonFullsimExtra, MuonMCTruthKinExtra, MuonPythiaExtra {};
struct MuonPowhegFullSimWTruth : MuonBase, MuonRecoExtra, MuonFullsimExtra, MuonMCTruthKinExtra, MuonPowhegExtra {};

struct MuonPythiaFullSimOverlayNoTruth : MuonBase, MuonRecoExtra, MuonFullsimExtra, MuonMCTruthKinExtra, MuonPythiaExtra, MuonPbPbExtra {};
struct MuonPowhegFullSimOverlayNoTruth : MuonBase, MuonRecoExtra, MuonFullsimExtra, MuonMCTruthKinExtra, MuonPowhegExtra, MuonPbPbExtra {};

struct MuonPythiaFullSimOverlayWTruth : MuonBase, MuonRecoExtra, MuonFullsimExtra, MuonMCTruthKinExtra, MuonPythiaExtra, MuonPbPbExtra {};
struct MuonPowhegFullSimOverlayWTruth : MuonBase, MuonRecoExtra, MuonFullsimExtra, MuonMCTruthKinExtra, MuonPowhegExtra, MuonPbPbExtra {};

