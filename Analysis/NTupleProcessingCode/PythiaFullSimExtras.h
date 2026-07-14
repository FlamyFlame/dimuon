#pragma once

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

    // Per-RECO-muon truth provenance + truth kinematics (bound only in store_mc_trigger mode).
    // These are the branches the reco-muon provenance classifier is built on (§3.0/D4):
    //   real     = prob > 0.5 && |muon_truth_id| == 13 && muon_truth_IsPrimary
    //   fake     = prob <= 0.5
    //   hadronic = |id| != 13 (punch-through) OR (|id| == 13 && !IsPrimary) (decay-in-flight)
    // The axis is (prob, |id|, IsPrimary) ONLY -- NEVER Pythia-signal-block membership, so a
    // real HIJING muon in the overlay is REAL (low_mass_dimuon_template_fit.md).
    std::vector<int>*    muon_truth_id        = nullptr;
    std::vector<bool>*   muon_truth_IsPrimary = nullptr;
    std::vector<float>*  muon_truth_pt        = nullptr;
    std::vector<float>*  muon_truth_eta       = nullptr;
    std::vector<float>*  muon_truth_phi       = nullptr;
    std::vector<int>*    muon_truth_charge    = nullptr;

    // provenance bookkeeping (store_mc_trigger; reported at end of run)
    long long n_prov_reco = 0, n_prov_real = 0, n_prov_fake = 0, n_prov_hadronic = 0;

    // trigger branches (bound only when store_mc_trigger; Run-3 chain names — the
    // trigger-enabled MC skims are Run-3 only). Mirrors DimuonDataAlgCoreT: per-muon
    // match = bare branch name (= mindR 0.02 nominal), pair-level 2mu4 = the
    // order-insensitive "_0_02" mindR branch, indexed by the skim's (i<j) pair block.
    std::vector<bool>*   muon_b_HLT_mu4        = nullptr; // muon_b_HLT_mu4_L1MU3V
    std::vector<bool>*   dimuon_b_2mu4_mindR   = nullptr; // dimuon_b_HLT_2mu4_L12MU3V_0_02
    std::vector<int>*    muon_pair_muon1_index = nullptr;
    std::vector<int>*    muon_pair_muon2_index = nullptr;
    Bool_t               b_HLT_mu4  {};                   // event-level decisions (diagnostics)
    Bool_t               b_HLT_2mu4 {};

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
    void ProcessEventFullsim(int ev_num);

    // store_mc_trigger ONLY (§3.0/D4): a RECO-SEEDED event loop over truth-matched REAL
    // muons, replacing the nominal truth-seeded Pythia-block loop. Kept as a separate
    // function so the nominal fullsim path (reco-efficiency, detector response,
    // template-fit MC) -- where a truth-seeded, Pythia-only denominator is the CORRECT
    // construction -- is byte-for-byte untouched.
    void ProcessEventFullsimMCTrig(int ev_num);

    // The real/hadronic/fake axis, on the RECO muon at index reco_ind.
    bool IsRealRecoMuon(int reco_ind);

    // The ONE definition of a reco muon's offline quantities + WP flags (shared by both loops).
    void FillRecoQuantities(muon_t& m, int reco_ind);
    void CheckBranchPtrsExtra();
    bool PassMuonMediumCuts(const muon_t& muon);

public:
    int  turn_on_track_charge = false;
    bool disable_ip_cut       = false; // debug: skip d0/z0 cuts in PassMuonMediumCuts
};
