template <class Derived>
class PowhegFullSimExtras{
protected:

    std::vector<float>   *muon_pair_muon1_pt           =nullptr;
    std::vector<float>   *muon_pair_muon1_eta       =nullptr;
    std::vector<float>   *muon_pair_muon1_phi          =nullptr;
    std::vector<int>     *muon_pair_muon1_ch         =nullptr;
    std::vector<int>     *muon_pair_muon1_bar         =nullptr;

    std::vector<float>   *muon_pair_muon2_pt       =nullptr;
    std::vector<float>   *muon_pair_muon2_eta          =nullptr;
    std::vector<float>   *muon_pair_muon2_phi         =nullptr;
    std::vector<int>     *muon_pair_muon2_ch         =nullptr;
    std::vector<int>     *muon_pair_muon2_bar         =nullptr;
	void InitInputExtra();
};