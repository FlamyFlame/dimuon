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
            return self().ResonanceTaggingImpl(!pair.same_sign, pair.minv, resonance_tagged_muon_index_list_reco);
        } else {
            std::cerr << "ResonanceTaggingReco requires PairT to have members `minv` and `same_sign`" << std::endl;
            return;
        }
    }
    void ResonanceTaggingTruth(){
        if (!self().mpairRef()) throw std::runtime_error("ResonanceTaggingTruth: mpair is nullptr!");
        auto& pair = *self().mpairRef();
        if constexpr (requires { pair.truth_minv; pair.truth_same_sign; }) {
            return self().ResonanceTaggingImpl(!pair.truth_same_sign, pair.truth_minv, resonance_tagged_muon_index_list_truth);
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
    void CheckBranchPtrsExtra();
    bool PassMuonMediumCuts(const muon_t& muon);

public:
    int turn_on_track_charge = false;
};
