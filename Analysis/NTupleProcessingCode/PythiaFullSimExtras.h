#pragma once

#include <vector>
#include "../MuonObjectsParamsAndHelpers/FullSimSampleType.h"

template <class PairT, class MuonT, class Derived>
class PythiaFullSimExtras {
    template <class, class, class, class...> friend class PythiaAlgCoreT;

protected:
    using pair_t = PairT;
    using muon_t = MuonT;

    float truth_match_prob_thrsh = 0.5;

    Derived& self() { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }

    // reco muon quantities
    std::vector<float>*  muon_pt           = nullptr;
    std::vector<float>*  muon_eta          = nullptr;
    std::vector<float>*  muon_phi          = nullptr;
    std::vector<int>*    muon_quality      = nullptr;
    std::vector<float>*  muon_deltaP_overP = nullptr;
    std::vector<float>*  muon_d0           = nullptr;
    std::vector<float>*  muon_z0           = nullptr;
    std::vector<float>*  muon_trk_pt       = nullptr;
    std::vector<float>*  muon_trk_eta      = nullptr;
    std::vector<float>*  muon_trk_phi      = nullptr;

    // reco-truth matching
    std::vector<float>*  muon_truth_prob    = nullptr;
    std::vector<int>*    muon_truth_barcode = nullptr;

    // trigger branches (bound only when store_mc_trigger; Run-3 chain names — the
    // trigger-enabled MC skims are Run-3 only). Mirrors DimuonDataAlgCoreT: per-muon
    // match = bare branch name (= mindR 0.02 nominal), pair-level 2mu4 = the
    // order-insensitive "_0_02" mindR branch, indexed by the skim's (i<j) pair block.
    std::vector<bool>*   muon_b_HLT_mu4        = nullptr; // muon_b_HLT_mu4_L1MU3V (full chain)
    std::vector<bool>*   muon_match_L1MU3V     = nullptr; // per-muon L1_MU3V RoI match (round-5 #3);
                                                          // only in RE-SKIMMED NTUPs, else nullptr
    bool                 has_l1_match = false;            // muon_match_L1MU3V present in the input
    std::vector<bool>*   dimuon_b_2mu4_mindR   = nullptr; // dimuon_b_HLT_2mu4_L12MU3V_0_02
    std::vector<int>*    muon_pair_muon1_index = nullptr;
    std::vector<int>*    muon_pair_muon2_index = nullptr;
    Bool_t               b_HLT_mu4  {};                   // event-level decisions (diagnostics)
    Bool_t               b_HLT_2mu4 {};

    // Reconstructed primary vertices (round 7, single-vertex sanity check). The skim dumps the
    // whole PrimaryVertices container unfiltered, so it ALWAYS contains one dummy beamspot
    // vertex with ntrk = 0 -- count only vtx_ntrk >= 2 (see Muon.h MuonFullsimExtra::n_vtx).
    // Bound only in store_mc_trigger mode; a skim without the branch leaves n_vtx = -1.
    std::vector<int>*    vtx_ntrk    = nullptr;
    bool                 has_vtx_ntrk = false;
    // Vertex z positions, in the same container order -- needed by the ALL-VERTEX
    // impact-parameter selection (Utilities/AllVertexIPSelection.h). Bound, together with
    // vtx_ntrk, on EVERY chain (not only in store_mc_trigger mode) whenever UseAllVertexIP().
    std::vector<float>*  vtx_z       = nullptr;

    // reco/ID efficiency SFs (skim tools; filled only for WP-passing muons, <=0 otherwise)
    std::vector<float>*  muon_eff_SF_medium    = nullptr;
    std::vector<float>*  muon_eff_SF_tight     = nullptr;
    // unfilled-SF bookkeeping: [0]=reco-matched muons, [1]=SF_medium<=0 (of which pass_medium),
    // [2]=SF_tight<=0 (of which pass_tight); split so the WP-passing fraction is exact
    long long n_sf_recomatched   = 0;
    long long n_sf_med_unfilled  = 0, n_sf_med_unfilled_wp  = 0, n_sf_med_wp  = 0;
    long long n_sf_tgt_unfilled  = 0, n_sf_tgt_unfilled_wp  = 0, n_sf_tgt_wp  = 0;
    void FinalizeExtra();

    // Look up the skim's pair index for two raw-NTUP reco-muon indices and return the
    // pair-level 2mu4 match decision. Throws if the pair is not in the skim block
    // (must not happen for two valid reco indices — fail fast rather than bias).
    bool LookupPairPass2mu4(int reco_ind_a, int reco_ind_b);

    // truth muon quantities (from truth_muon_* branches)
    std::vector<float>*  truth_muon_pt  = nullptr;
    std::vector<float>*  truth_muon_eta = nullptr;
    std::vector<float>*  truth_muon_phi = nullptr;
    std::vector<int>*    truth_muon_ch  = nullptr;

    std::vector<int> resonance_tagged_muon_index_list_reco {};
    std::vector<int> resonance_tagged_muon_index_list_truth {};

    void ResonanceTaggingReco(){
        if (!self().mpairRef()) throw std::runtime_error("ResonanceTaggingReco: mpair is nullptr!");
        auto& pair = *self().mpairRef();
        if constexpr (requires { pair.minv; pair.same_sign; }) {
            return self().ResonanceTaggingImpl(!pair.same_sign, pair.minv, resonance_tagged_muon_index_list_reco, self().pmsRef().minv_cuts_v2);
        } else {
            std::cerr << "ResonanceTaggingReco requires PairT to have members `minv` and `same_sign`" << std::endl;
            return;
        }
    }
    void ResonanceTaggingTruth(){
        if (!self().mpairRef()) throw std::runtime_error("ResonanceTaggingTruth: mpair is nullptr!");
        auto& pair = *self().mpairRef();
        if constexpr (requires { pair.truth_minv; pair.truth_same_sign; }) {
            return self().ResonanceTaggingImpl(!pair.truth_same_sign, pair.truth_minv, resonance_tagged_muon_index_list_truth, self().pmsRef().minv_cuts_v2);
        } else {
            std::cerr << "ResonanceTaggingTruth requires PairT to have members `truth_minv` and `truth_same_sign`" << std::endl;
            return;
        }
    }

    void InitInputExtra();
    void InitParamsExtra(){
        self().setIsFullsim(true);
    }
    // Truth-seeded, Pythia-ONLY event loop. Nominal fullsim (reco-eff, det-response,
    // template-fit MC) AND store_mc_trigger (MC trigger efficiency) BOTH use it: only
    // Pythia-truth muons matched to a reco muon may enter, so in the HIJING overlay the
    // real HIJING muons are excluded (they carry a different event weight than Pythia and
    // belong only to the template-fit background). store_mc_trigger differs only by the
    // per-muon/per-pair trigger fields, the RECO-loose single-muon gate, and no truth-pT gate
    // on the singles turn-on. (mc_trigger_efficiency.md round-5 change #1: this REVERTS the
    // 2026-07-14 reco-seeded ProcessEventFullsimMCTrig / IsRealRecoMuon, now deleted.)
    void ProcessEventFullsim(int ev_num);

    // The ONE definition of a reco muon's offline quantities + WP flags.
    void FillRecoQuantities(muon_t& m, int reco_ind);
    void CheckBranchPtrsExtra();
    bool PassMuonMediumCuts(const muon_t& muon);

    // ---- ALL-VERTEX impact-parameter selection (pp-CONDITIONS fullsim ONLY) -------------
    // Mirrors the pp DATA selection (NTupleProcessingCode/PPExtras.c) through the shared
    // Utilities/AllVertexIPSelection.h, so eps_reco is measured on the SAME pair selection the
    // data applies. Without this the pp24 fullsim -- which has genuine pile-up, track-bearing
    // vertices 1:5.0%, 2:14.7%, 3:21.8%, 4:22.7%, essentially the same as data's 1:6.4%,
    // 2:17.9%, 3:25.1%, 4:23.0% -- would supply an efficiency measured on a ~1.6-1.9 % narrower
    // selection than the data it corrects.
    // ONLY FullSimSampleType::pp. The HIJING overlay has exactly one track-bearing vertex in
    // 100 % of events, so it is an exact no-op there and is left off explicitly rather than
    // relied on; zmumu / noovl / data likewise keep the primary-vertex cut.
    bool UseAllVertexIP() const {
        return self().getFullSimSampleType() == FullSimSampleType::pp;
    }
    int  PairAllVertexIndex(const muon_t& m1, const muon_t& m2, bool* pass_primary_out);

    // Counters for the end-of-job report (no output branch needed).
    long long n_allvtx_pairs_pass{0}, n_allvtx_pairs_secondary{0}, n_allvtx_pairs_added{0};

public:
    int  turn_on_track_charge = false;
    bool disable_ip_cut       = false; // debug: skip d0/z0 cuts in PassMuonMediumCuts
};
