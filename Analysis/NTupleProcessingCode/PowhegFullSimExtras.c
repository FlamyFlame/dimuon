#include "PowhegFullSimExtras.h"
#include "../Utilities/tchain_helpers.h"

template <class PairT, class MuonT, class Derived>
void PowhegFullSimExtras<PairT, MuonT, Derived>::InitInputExtra(){
    enable_and_bind(self().fChainRef(), "muon_pt"                , &muon_pt);
    enable_and_bind(self().fChainRef(), "muon_eta"               , &muon_eta);
    enable_and_bind(self().fChainRef(), "muon_phi"               , &muon_phi);
    enable_and_bind(self().fChainRef(), "muon_quality"           , &muon_quality);
    enable_and_bind(self().fChainRef(), "muon_deltaP_overP"      , &muon_deltaP_overP);
    enable_and_bind(self().fChainRef(), "muon_d0"                , &muon_d0);
    enable_and_bind(self().fChainRef(), "muon_z0"                , &muon_z0);
    enable_and_bind(self().fChainRef(), "muon_trk_pt"            , &muon_trk_pt);
    enable_and_bind(self().fChainRef(), "muon_trk_eta"           , &muon_trk_eta);
    enable_and_bind(self().fChainRef(), "muon_trk_phi"           , &muon_trk_phi);
    enable_and_bind(self().fChainRef(), "muon_truth_prob"        , &muon_truth_prob);
    enable_and_bind(self().fChainRef(), "muon_truth_barcode"     , &muon_truth_barcode);
    enable_and_bind(self().fChainRef(), "truth_muon_pt"          , &truth_muon_pt);
    enable_and_bind(self().fChainRef(), "truth_muon_eta"         , &truth_muon_eta);
    enable_and_bind(self().fChainRef(), "truth_muon_phi"         , &truth_muon_phi);
    enable_and_bind(self().fChainRef(), "truth_muon_ch"          , &truth_muon_ch);
}

template <class PairT, class MuonT, class Derived>
bool PowhegFullSimExtras<PairT, MuonT, Derived>::PassMuonMediumCuts(const muon_t& muon){
    // apply all cuts applied to data muons --> obtain one reco efficiency for analysis

    // quality cuts
    if ((muon.quality&1  )==0) return false;//combined muon
    if ((muon.quality&8  )==0) return false;//Medium muon
    if ((muon.quality&32 )==0) return false;//IDCuts
    if ((muon.quality&256)==0) return false;//MuonCuts

    // reco kinematic cuts
    if ((muon.eta) > 2.4) return false;
    if (muon.pt < 4) return false;

    // HF muon cut    
    if (muon.dP_overP > self().pmsRef().deltaP_overP_thrsh ) return false;
    
    //cut on d0 & z0 sin(theta) against fake muons
    double z0sinTheta = fabs(muon.z0 * sin(2.0*atan(exp(-muon.eta))));
    bool pass_d0_z0_cuts = (fabs(muon.d0) < self().pmsRef().d0cut && z0sinTheta < self().pmsRef().z0cut);
    if (!pass_d0_z0_cuts) return false;

    // if track charge saved, require muon + track charge to agree
    if (turn_on_track_charge){
        if (muon.trk_charge != muon.charge) return false;
    }
    
    return true;
}

template <class PairT, class MuonT, class Derived>
void PowhegFullSimExtras<PairT, MuonT, Derived>::CheckBranchPtrsExtra(){
    if (self().debug_mode) std::cout << "Calling PowhegFullSimExtras::CheckBranchPtrsExtra" << std::endl;
    
    auto require = [&](auto* p, const char* name){
        if(!p) throw std::runtime_error(std::string("Null branch pointer: ") + name);
    };

    require(muon_pt, "muon_pt");
    require(muon_eta, "muon_eta");
    require(muon_phi, "muon_phi");
    require(muon_quality, "muon_quality");
    require(muon_deltaP_overP, "muon_deltaP_overP");
    require(muon_d0, "muon_d0");
    require(muon_z0, "muon_z0");
    require(muon_trk_pt, "muon_trk_pt");
    require(muon_trk_eta, "muon_trk_eta");
    require(muon_trk_phi, "muon_trk_phi");
    require(muon_truth_prob, "muon_truth_prob");
    require(muon_truth_barcode, "muon_truth_barcode");
    require(truth_muon_pt, "truth_muon_pt");
    require(truth_muon_eta, "truth_muon_eta");
    require(truth_muon_phi, "truth_muon_phi");
    require(truth_muon_ch, "truth_muon_ch");

    if (self().debug_mode) std::cout << "muon_pt size " << muon_pt->size() << ", " << muon_pt->at(0) << std::endl;
    if (self().debug_mode) std::cout << "truth_muon_eta size " << truth_muon_eta->size() << ", " << truth_muon_eta->at(0) << std::endl;
}

template <class PairT, class MuonT, class Derived>
void PowhegFullSimExtras<PairT, MuonT, Derived>::ProcessEventFullsim(int ev_num){ // per-event analysis    
    
    // -------- save barcode for truth muons that match with reco muons with prob > truth_match_prob_thrsh --------
    std::vector<int> real_muon_truth_barcode_list({}); // truth barcodes of "real" reco muons matched to truth muons with probability > threshold
    std::vector<int> fake_muon_ind_list({});

    if (muon_truth_barcode->size() != muon_truth_prob->size()){
        throw std::runtime_error("muon_truth_barcode & muon_truth_prob must have the same size!");
    }

    // loop over reco muons + see if they are matched to truth muons with probability > threshold
    for (int ind = 0; ind < muon_truth_barcode->size(); ind++){
        if (muon_truth_prob->at(ind) > truth_match_prob_thrsh){ // reco muon matched to truth muon --> real muon
            real_muon_truth_barcode_list.push_back(muon_truth_barcode->at(ind));
        } else { // fake muon
            fake_muon_ind_list.push_back(ind);
        }
    }

    if (muon_truth_barcode->size() != muon_pt->size()){
        throw std::runtime_error("Muon truth (including fakes) & muon vectors (muon_truth_barcode & muon_pt) must have the same size!");
    }

    float event_crossx = (self().EventWeightsRef()->size() == 0)? 1. : self().EventWeightsRef()->at(0);
    float event_weight = (self().EventWeightsRef()->size() == 0)? 1. : event_crossx * self().filter_effcy;

    std::vector<muon_t> truth_muon_list ({});

    for (int truth_ind = 0; truth_ind < truth_muon_pt->size(); truth_ind++){
        // -------- truth quantities --------
        muon_t cur_muon;
        cur_muon.ind = truth_ind;
        cur_muon.ev_num = ev_num;
        cur_muon.ev_weight = event_weight;

        cur_muon.truth_pt       = fabs(truth_muon_pt->at(truth_ind))/1000.0;
        cur_muon.truth_eta      = truth_muon_eta->at(truth_ind);
        cur_muon.truth_phi      = truth_muon_phi->at(truth_ind);
        cur_muon.truth_charge   = truth_muon_ch->at(truth_ind);
        cur_muon.truth_bar      = self().truth_muon_barcode->at(truth_ind);

        // -------- truth-to-reco-muon matching --------
        int reco_ind = -1;
        auto it = std::find(real_muon_truth_barcode_list.begin(), real_muon_truth_barcode_list.end(), cur_muon.truth_bar);

        if (it != real_muon_truth_barcode_list.end()) { // found reco muon
            cur_muon.reco_match = true;
            reco_ind = std::distance(real_muon_truth_barcode_list.begin(), it);

            // fill reco quantities
            cur_muon.pt = fabs(muon_pt->at(reco_ind))/1000.0;
            cur_muon.eta = muon_eta->at(reco_ind);
            cur_muon.phi = muon_phi->at(reco_ind);
            cur_muon.charge = (muon_pt ->at(reco_ind) > 0)? 1:-1;
            
            cur_muon.dP_overP = muon_deltaP_overP->at(reco_ind);
            cur_muon.z0 = muon_z0->at(reco_ind);
            cur_muon.d0 = muon_d0->at(reco_ind);
            cur_muon.quality = muon_quality->at(reco_ind);
            
            cur_muon.trk_pt      = fabs(muon_trk_pt->at(reco_ind))/1000.0;
            cur_muon.trk_eta     = muon_trk_eta->at(reco_ind);
            cur_muon.trk_phi     = muon_trk_phi->at(reco_ind);

            if (turn_on_track_charge){
                cur_muon.trk_charge = (muon_trk_pt ->at(reco_ind) > 0)? 1:-1;//sign of pt stores charge
            } else{
                cur_muon.trk_charge = 0;
            }

            cur_muon.pass_medium = PassMuonMediumCuts(cur_muon);
            cur_muon.pass_tight = (cur_muon.pass_medium && (cur_muon.quality&16));

        } else { // no reco muon found
            cur_muon.reco_match = false;
            cur_muon.pass_medium = false;
            cur_muon.pass_tight = false;
        }

        if(self().output_single_muon_tree){
            self().muon_raw_ptr = &cur_muon;
            self().FillSingleMuonTree();
        }
        
        truth_muon_list.push_back(std::move(cur_muon));
    } // end loop over single truth muons

    if(self().output_single_muon_tree) return;

    // -------- build truth pairs --------

    self().muon_pair_list_cur_event_pre_resonance_cut.clear();
    self().resonance_tagged_muon_index_list_reco.clear(); // MUST CLEAR for each event!!
    self().resonance_tagged_muon_index_list_truth.clear(); // MUST CLEAR for each event!!

    for (int i = 0; i < truth_muon_list.size() - 1; i++){ // begin first loop over muon pairs
        for (int j = i+1; j < truth_muon_list.size(); j++){
            if (self().mpairRef()) self().mpairRef()->Clear();
            else self().mpairRef() = std::make_shared<pair_t>();

            self().mpairRef()->weight     = event_weight;
            self().mpairRef()->crossx     = event_crossx;

            self().mpairRef()->m1 = truth_muon_list.at(i); // default copy constructor
            self().mpairRef()->m2 = truth_muon_list.at(j); // default copy constructor

            self().h_cutAcceptanceRef()[self().mpairRef()->m1.truth_charge != self().mpairRef()->m2.truth_charge]->Fill(double(nocut) + 0.5, self().mpairRef()->weight);
    
            if (abs(self().mpairRef()->weight) > self().crossx_cut * self().filter_effcy) continue;

            if (!self().PassCuts()) continue; // truth cuts

            self().mpairRef()->Update(); // only afterwards can we use mpairRef()->same_sign
            
            if (self().mpairRef()->m1.reco_match && self().mpairRef()->m2.reco_match){                
                self().ResonanceTaggingReco();
                self().ResonanceTaggingTruth();
                self().mpairRef()->pair_pass_medium = (self().mpairRef()->m1.pass_medium && self().mpairRef()->m2.pass_medium);
                self().mpairRef()->pair_pass_tight  = (self().mpairRef()->m1.pass_tight  && self().mpairRef()->m2.pass_tight);
            }else{
                self().mpairRef()->pair_pass_medium = false;
                self().mpairRef()->pair_pass_tight  = false;
            }

            self().muon_pair_list_cur_event_pre_resonance_cut.push_back(std::move(self().mpairRef()));
        }
    } // end 1st loop over truth muon pairs

    for(int pair_ind = 0; pair_ind < self().muon_pair_list_cur_event_pre_resonance_cut.size(); pair_ind++){//second loop over all muon-pairs in the event
        self().mpairRef() = std::move(self().muon_pair_list_cur_event_pre_resonance_cut.at(pair_ind));

        if (!self().mpairRef()){
            std::cerr << "mpairRef() at second muon-pair loop NOT found! Skipping pair in output tree filling." << std::endl;
            continue;
        }

        std::vector<int>::iterator itres_m1;
        std::vector<int>::iterator itres_m2;

        if (self().mpairRef()->m1.reco_match && self().mpairRef()->m2.reco_match){            
            self().mpairRef()->pair_pass_resonance_reco = true;
            self().mpairRef()->pair_pass_resonance_truth = true;

            itres_m1 = std::find(self().resonance_tagged_muon_index_list_reco.begin(),self().resonance_tagged_muon_index_list_reco.end(),self().mpairRef()->m1.ind);
            self().mpairRef()->pair_pass_resonance_reco &= (itres_m1 == self().resonance_tagged_muon_index_list_reco.end());

            itres_m2 = std::find(self().resonance_tagged_muon_index_list_reco.begin(),self().resonance_tagged_muon_index_list_reco.end(),self().mpairRef()->m2.ind);
            self().mpairRef()->pair_pass_resonance_reco &= (itres_m2 == self().resonance_tagged_muon_index_list_reco.end());

            itres_m1 = std::find(resonance_tagged_muon_index_list_truth.begin(),resonance_tagged_muon_index_list_truth.end(),self().mpairRef()->m1.ind);
            self().mpairRef()->pair_pass_resonance_truth &= (itres_m1 == resonance_tagged_muon_index_list_truth.end());

            itres_m2 = std::find(resonance_tagged_muon_index_list_truth.begin(),resonance_tagged_muon_index_list_truth.end(),self().mpairRef()->m2.ind);
            self().mpairRef()->pair_pass_resonance_truth &= (itres_m2 == resonance_tagged_muon_index_list_truth.end());

            self().mpairRef()->pair_pass_medium_and_resonance = self().mpairRef()->pair_pass_resonance_reco && self().mpairRef()->pair_pass_medium;
            self().mpairRef()->pair_pass_tight_and_resonance  = self().mpairRef()->pair_pass_resonance_reco && self().mpairRef()->pair_pass_tight;
        }else{
            self().mpairRef()->pair_pass_resonance_reco = false;
            self().mpairRef()->pair_pass_resonance_truth = false;
            self().mpairRef()->pair_pass_medium = false;
            self().mpairRef()->pair_pass_tight  = false;
            self().mpairRef()->pair_pass_medium_and_resonance = false;
            self().mpairRef()->pair_pass_tight_and_resonance  = false;
        }

        if (self().getPerformTruth()) self().PerformTruthPairAnalysisHook();

        self().FillMuonPairTree();
    } // end 2nd loop over truth muon pairs
}
