template <class PairT, class MuonT, class Derived>
class PowhegFullSimExtras {
    template <class, class, class, class...> friend class PowhegAlgCoreT;

protected:
    using pair_t = PairT;
    using muon_t = MuonT;

    float truth_match_prob_thrsh = 0.5;

    Derived& self() { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }

    // reco muon quantities
    std::vector<float>*      muon_pt;
    std::vector<float>*      muon_eta;
    std::vector<float>*      muon_phi;
    std::vector<int>*        muon_quality;
    std::vector<float>*      muon_deltaP_overP;
    std::vector<float>*      muon_d0;
    std::vector<float>*      muon_z0;
    std::vector<float>*      muon_trk_pt;
    std::vector<float>*      muon_trk_eta;
    std::vector<float>*      muon_trk_phi;

    // reco-truth matching
    std::vector<float>*      muon_truth_prob;
    std::vector<int>*        muon_truth_barcode;

    // truth muon quantities
    std::vector<float>*      truth_muon_pt;
    std::vector<float>*      truth_muon_eta;
    std::vector<float>*      truth_muon_phi;
    std::vector<float>*      truth_muon_ch;
    std::vector<float>*      truth_muon_barcode;

	void InitInputExtra();
    void InitParamsExtra(){
        self().setIsFullsim(true);
    }
    void ProcessEventFullsim(int ev_num);
    bool PassMuonMediumCuts(const muon_t& muon);

public:
    int turn_on_track_charge = false;
};