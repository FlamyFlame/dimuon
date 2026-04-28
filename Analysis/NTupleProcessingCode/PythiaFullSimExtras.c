#include "PythiaFullSimExtras.h"
#include "../Utilities/tchain_helpers.h"

template <class PairT, class MuonT, class Derived>
void PythiaFullSimExtras<PairT, MuonT, Derived>::InitInputExtra(){
    // Bind reco branches on every chain in evChains_kn_beam
    auto bind_reco = [&](TChain* ch) {
        if (!ch) return;
        enable_and_bind(ch, "muon_pt"           , &muon_pt);
        enable_and_bind(ch, "muon_eta"          , &muon_eta);
        enable_and_bind(ch, "muon_phi"          , &muon_phi);
        enable_and_bind(ch, "muon_quality"      , &muon_quality);
        enable_and_bind(ch, "muon_deltaP_overP" , &muon_deltaP_overP);
        enable_and_bind(ch, "muon_d0"           , &muon_d0);
        enable_and_bind(ch, "muon_z0"           , &muon_z0);
        enable_and_bind(ch, "muon_trk_pt"       , &muon_trk_pt);
        enable_and_bind(ch, "muon_trk_eta"      , &muon_trk_eta);
        enable_and_bind(ch, "muon_trk_phi"      , &muon_trk_phi);
        enable_and_bind(ch, "muon_truth_prob"   , &muon_truth_prob);
        enable_and_bind(ch, "muon_truth_barcode", &muon_truth_barcode);
        enable_and_bind(ch, "truth_muon_pt"     , &truth_muon_pt);
        enable_and_bind(ch, "truth_muon_eta"    , &truth_muon_eta);
        enable_and_bind(ch, "truth_muon_phi"    , &truth_muon_phi);
        enable_and_bind(ch, "truth_muon_ch"     , &truth_muon_ch);
    };

    for (int ikin = 0; ikin < self().nKinRanges; ikin++)
        for (int ibeam = 0; ibeam < self().nBeamTypes; ibeam++)
            bind_reco(self().evChains_kn_beam[ikin][ibeam]);
}

template <class PairT, class MuonT, class Derived>
bool PythiaFullSimExtras<PairT, MuonT, Derived>::PassMuonMediumCuts(const muon_t& muon){
    if ((muon.quality&1  )==0) return false; // combined muon
    if ((muon.quality&8  )==0) return false; // Medium muon
    if ((muon.quality&32 )==0) return false; // IDCuts
    if ((muon.quality&256)==0) return false; // MuonCuts

    if (fabs(muon.eta) > 2.4) return false;
    if (muon.pt < 4) return false;

    if (muon.dP_overP > self().pmsRef().deltaP_overP_thrsh) return false;

    if (!self().disable_ip_cut) {
        double z0sinTheta = fabs(muon.z0 * sin(2.0*atan(exp(-muon.eta))));
        if (fabs(muon.d0) >= self().pmsRef().d0cut || z0sinTheta >= self().pmsRef().z0cut) return false;
    }

    if (turn_on_track_charge){
        if (muon.trk_charge != muon.charge) return false;
    }
    return true;
}

template <class PairT, class MuonT, class Derived>
void PythiaFullSimExtras<PairT, MuonT, Derived>::CheckBranchPtrsExtra(){
    if (self().debug_mode) std::cout << "Calling PythiaFullSimExtras::CheckBranchPtrsExtra" << std::endl;

    auto require = [&](auto* p, const char* name){
        if (!p) throw std::runtime_error(std::string("Null branch pointer: ") + name);
    };

    require(muon_pt,            "muon_pt");
    require(muon_eta,           "muon_eta");
    require(muon_phi,           "muon_phi");
    require(muon_quality,       "muon_quality");
    require(muon_deltaP_overP,  "muon_deltaP_overP");
    require(muon_d0,            "muon_d0");
    require(muon_z0,            "muon_z0");
    require(muon_trk_pt,        "muon_trk_pt");
    require(muon_trk_eta,       "muon_trk_eta");
    require(muon_trk_phi,       "muon_trk_phi");
    require(muon_truth_prob,    "muon_truth_prob");
    require(muon_truth_barcode, "muon_truth_barcode");
    require(truth_muon_pt,      "truth_muon_pt");
    require(truth_muon_eta,     "truth_muon_eta");
    require(truth_muon_phi,     "truth_muon_phi");
    require(truth_muon_ch,      "truth_muon_ch");
    require(self().truth_muon_barcode, "truth_muon_barcode (core)");
}

template <class PairT, class MuonT, class Derived>
void PythiaFullSimExtras<PairT, MuonT, Derived>::ProcessEventFullsim(int ev_num){

    // ---- Build list of "real" reco muon truth barcodes (prob > threshold) ----
    std::vector<int> real_muon_truth_barcode_list;
    std::vector<int> fake_muon_ind_list;

    if (muon_truth_barcode->size() != muon_truth_prob->size())
        throw std::runtime_error("muon_truth_barcode & muon_truth_prob size mismatch!");

    for (int ind = 0; ind < (int)muon_truth_barcode->size(); ind++){
        if (muon_truth_prob->at(ind) > truth_match_prob_thrsh)
            real_muon_truth_barcode_list.push_back(muon_truth_barcode->at(ind));
        else
            fake_muon_ind_list.push_back(ind);
    }

    if (muon_truth_barcode->size() != muon_pt->size())
        throw std::runtime_error("muon_truth_barcode & muon_pt size mismatch!");

    // Event weight: precomputed per beam/kn as ami_weight * isospin_ratio / N_beam
    float event_weight = static_cast<float>(self().fullsim_weight_factor);

    // ---- Build truth muon list with optional reco matching ----
    std::vector<muon_t> truth_muon_list;

    for (int truth_ind = 0; truth_ind < (int)truth_muon_pt->size(); truth_ind++){
        muon_t cur_muon;
        cur_muon.ind      = truth_ind;
        cur_muon.ev_num   = ev_num;
        cur_muon.ev_weight = event_weight;

        cur_muon.truth_pt     = fabs(truth_muon_pt->at(truth_ind)) / 1000.0f;
        cur_muon.truth_eta    = truth_muon_eta->at(truth_ind);
        cur_muon.truth_phi    = truth_muon_phi->at(truth_ind);
        cur_muon.truth_charge = truth_muon_ch->at(truth_ind);
        cur_muon.truth_bar    = self().truth_muon_barcode->at(truth_ind);

        // Truth-to-reco matching via barcode
        auto it = std::find(real_muon_truth_barcode_list.begin(),
                            real_muon_truth_barcode_list.end(),
                            cur_muon.truth_bar);

        if (it != real_muon_truth_barcode_list.end()){
            cur_muon.reco_match = true;
            int reco_ind = std::distance(real_muon_truth_barcode_list.begin(), it);

            cur_muon.pt      = fabs(muon_pt->at(reco_ind)) / 1000.0f;
            cur_muon.eta     = muon_eta->at(reco_ind);
            cur_muon.phi     = muon_phi->at(reco_ind);
            cur_muon.charge  = (muon_pt->at(reco_ind) > 0) ? 1 : -1;

            cur_muon.dP_overP = muon_deltaP_overP->at(reco_ind);
            cur_muon.z0       = muon_z0->at(reco_ind);
            cur_muon.d0       = muon_d0->at(reco_ind);
            cur_muon.quality  = muon_quality->at(reco_ind);

            cur_muon.trk_pt  = fabs(muon_trk_pt->at(reco_ind)) / 1000.0f;
            cur_muon.trk_eta = muon_trk_eta->at(reco_ind);
            cur_muon.trk_phi = muon_trk_phi->at(reco_ind);

            if (turn_on_track_charge)
                cur_muon.trk_charge = (muon_trk_pt->at(reco_ind) > 0) ? 1 : -1;
            else
                cur_muon.trk_charge = 0;

            cur_muon.pass_medium = PassMuonMediumCuts(cur_muon);
            cur_muon.pass_tight  = (cur_muon.pass_medium && (cur_muon.quality & 16));
        } else {
            cur_muon.reco_match  = false;
            cur_muon.pass_medium = false;
            cur_muon.pass_tight  = false;
        }

        if constexpr (requires { self().FillMuonOverlay(cur_muon); }) {
            self().FillMuonOverlay(cur_muon);
        }

        if (self().output_single_muon_tree){
            if (cur_muon.truth_pt > 4.0 && fabs(cur_muon.truth_eta) < 2.4){
                self().muon_raw_ptr = &cur_muon;
                self().FillSingleMuonTree();
            }
        }

        truth_muon_list.push_back(std::move(cur_muon));
    }

    if (self().output_single_muon_tree) return;

    // ---- Build truth pairs ----
    self().muon_pair_list_cur_event_pre_resonance_cut.clear();
    self().resonance_tagged_muon_index_list_reco.clear();
    self().resonance_tagged_muon_index_list_truth.clear();

    for (int i = 0; i < (int)truth_muon_list.size() - 1; i++){
        for (int j = i + 1; j < (int)truth_muon_list.size(); j++){
            if (self().mpairRef()) self().mpairRef()->Clear();
            else self().mpairRef() = std::make_shared<pair_t>();

            self().mpairRef()->weight  = event_weight;
            self().mpairRef()->crossx  = event_weight;

            self().mpairRef()->m1 = truth_muon_list.at(i);
            self().mpairRef()->m2 = truth_muon_list.at(j);

            self().h_cutAcceptanceRef()[self().mpairRef()->m1.truth_charge != self().mpairRef()->m2.truth_charge]
                ->Fill(double(nocut) + 0.5, event_weight);

            if (!self().PassCuts()) continue; // truth fiducial cuts

            self().mpairRef()->Update(); // compute reco+truth kinematics

            if constexpr (requires { self().FillPairOverlay(); }) {
                self().FillPairOverlay();
            }

            if (self().mpairRef()->m1.reco_match && self().mpairRef()->m2.reco_match){
                ResonanceTaggingReco();
                ResonanceTaggingTruth();
                self().mpairRef()->pair_pass_medium = (self().mpairRef()->m1.pass_medium && self().mpairRef()->m2.pass_medium);
                self().mpairRef()->pair_pass_tight  = (self().mpairRef()->m1.pass_tight  && self().mpairRef()->m2.pass_tight);
            } else {
                self().mpairRef()->pair_pass_medium = false;
                self().mpairRef()->pair_pass_tight  = false;
            }

            self().muon_pair_list_cur_event_pre_resonance_cut.push_back(std::move(self().mpairRef()));
        }
    }

    // ---- Second loop: apply resonance veto and fill output ----
    for (int pair_ind = 0; pair_ind < (int)self().muon_pair_list_cur_event_pre_resonance_cut.size(); pair_ind++){
        self().mpairRef() = std::move(self().muon_pair_list_cur_event_pre_resonance_cut.at(pair_ind));
        if (!self().mpairRef()){
            std::cerr << "mpairRef() at second muon-pair loop NOT found!" << std::endl;
            continue;
        }

        if (self().mpairRef()->m1.reco_match && self().mpairRef()->m2.reco_match){
            self().mpairRef()->pair_pass_resonance_reco  = true;
            self().mpairRef()->pair_pass_resonance_truth = true;

            auto it1 = std::find(resonance_tagged_muon_index_list_reco.begin(),
                                 resonance_tagged_muon_index_list_reco.end(),
                                 self().mpairRef()->m1.ind);
            self().mpairRef()->pair_pass_resonance_reco &= (it1 == resonance_tagged_muon_index_list_reco.end());

            auto it2 = std::find(resonance_tagged_muon_index_list_reco.begin(),
                                 resonance_tagged_muon_index_list_reco.end(),
                                 self().mpairRef()->m2.ind);
            self().mpairRef()->pair_pass_resonance_reco &= (it2 == resonance_tagged_muon_index_list_reco.end());

            it1 = std::find(resonance_tagged_muon_index_list_truth.begin(),
                            resonance_tagged_muon_index_list_truth.end(),
                            self().mpairRef()->m1.ind);
            self().mpairRef()->pair_pass_resonance_truth &= (it1 == resonance_tagged_muon_index_list_truth.end());

            it2 = std::find(resonance_tagged_muon_index_list_truth.begin(),
                            resonance_tagged_muon_index_list_truth.end(),
                            self().mpairRef()->m2.ind);
            self().mpairRef()->pair_pass_resonance_truth &= (it2 == resonance_tagged_muon_index_list_truth.end());

            self().mpairRef()->pair_pass_medium_and_resonance =
                self().mpairRef()->pair_pass_resonance_reco && self().mpairRef()->pair_pass_medium;
            self().mpairRef()->pair_pass_tight_and_resonance =
                self().mpairRef()->pair_pass_resonance_reco && self().mpairRef()->pair_pass_tight;
        } else {
            self().mpairRef()->pair_pass_resonance_reco  = false;
            self().mpairRef()->pair_pass_resonance_truth = false;
            self().mpairRef()->pair_pass_medium                = false;
            self().mpairRef()->pair_pass_tight                 = false;
            self().mpairRef()->pair_pass_medium_and_resonance  = false;
            self().mpairRef()->pair_pass_tight_and_resonance   = false;
        }

        if (self().getPerformTruth()) self().PerformTruthPairAnalysisHook();

        self().FillMuonPairTree();
    }
}
