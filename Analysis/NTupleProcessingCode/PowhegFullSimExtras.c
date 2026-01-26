#include "PowhegFullSimExtras.h"

template <class Derived>
void PowhegFullSimExtras<Derived>::InitInputExtra(){
    fChain()->SetBranchAddress("muon_pt"                , &muon_pt);
    fChain()->SetBranchAddress("muon_eta"               , &muon_eta);
    fChain()->SetBranchAddress("muon_phi"               , &muon_phi);
    fChain()->SetBranchAddress("muon_quality"           , &muon_quality);
    fChain()->SetBranchAddress("muon_deltaP_overP"      , &muon_deltaP_overP);
    fChain()->SetBranchAddress("muon_d0"                , &muon_d0);
    fChain()->SetBranchAddress("muon_z0"                , &muon_z0);
    fChain()->SetBranchAddress("muon_trk_pt"            , &muon_trk_pt);
    fChain()->SetBranchAddress("muon_trk_eta"           , &muon_trk_eta);
    fChain()->SetBranchAddress("muon_trk_phi"           , &muon_trk_phi);
    fChain()->SetBranchAddress("muon_truth_prob"        , &muon_truth_prob);
    fChain()->SetBranchAddress("muon_truth_barcode"     , &muon_truth_barcode);
    fChain()->SetBranchAddress("truth_muon_pt"          , &truth_muon_pt);
    fChain()->SetBranchAddress("truth_muon_eta"         , &truth_muon_eta);
    fChain()->SetBranchAddress("truth_muon_phi"         , &truth_muon_phi);
    fChain()->SetBranchAddress("truth_muon_ch"          , &truth_muon_ch);
    fChain()->SetBranchAddress("truth_muon_barcode"     , &truth_muon_barcode);

    fChain()->SetBranchStatus("muon_pt"                 ,1);
    fChain()->SetBranchStatus("muon_eta"                ,1);
    fChain()->SetBranchStatus("muon_phi"                ,1);
    fChain()->SetBranchStatus("muon_quality"            ,1);
    fChain()->SetBranchStatus("muon_deltaP_overP"       ,1);
    fChain()->SetBranchStatus("muon_d0"                 ,1);
    fChain()->SetBranchStatus("muon_z0"                 ,1);
    fChain()->SetBranchStatus("muon_trk_pt"             ,1);
    fChain()->SetBranchStatus("muon_trk_eta"            ,1);
    fChain()->SetBranchStatus("muon_trk_phi"            ,1);
    fChain()->SetBranchStatus("muon_truth_prob"         ,1);
    fChain()->SetBranchStatus("muon_truth_barcode"      ,1);
    fChain()->SetBranchStatus("truth_muon_pt"           ,1);
    fChain()->SetBranchStatus("truth_muon_eta"          ,1);
    fChain()->SetBranchStatus("truth_muon_phi"          ,1);
    fChain()->SetBranchStatus("truth_muon_ch"           ,1);
    fChain()->SetBranchStatus("truth_muon_barcode"      ,1);
}

template <class PairT, class MuonT, class Derived, class... Extras>
bool PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::PassMuonMediumCuts(const muon_t& muon){
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
    if (muon.dP_overP > self().pms.deltaP_overP_thrsh ) return false;
    
    //cut on d0 & z0 sin(theta) against fake muons
    double z0sinTheta = fabs(muon.z0 * sin(2.0*atan(exp(-muon.eta))));
    bool pass_d0_z0_cuts = (fabs(muon.d0) < self().pms.d0cut && z0sinTheta < self().pms.z0cut);
    if (!pass_d0_z0_cuts) return false;

    // if track charge saved, require muon + track charge to agree
    if (turn_on_track_charge){
        if (muon.trk_charge != muon.charge) return false;
    }
    
    return true;
}


template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::ProcessEventFullsim(int ev_num){ // per-event analysis
    
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
            fake_muon_ind_list.push_back()
        }
    }

    if (muon_truth_barcode->size() != muon_pt->size()){
        throw std::runtime_error("Muon truth (including nulls) & muon vectors (muon_truth_barcode & muon_pt) must have the same size!");
    }

    std::vector<muon_t> truth_muon_list ({});

    for (int muon_ind = 0; muon_ind < truth_muon_pt->size(); muon_ind++){
        // -------- truth quantities --------
        muon_t cur_muon;
        cur_muon.ind = muon_ind;
        cur_muon.ev_num = ev_num;

        cur_muon.truth_pt       = fabs(truth_muon_pt->at(muon_ind))/1000.0;
        cur_muon.truth_eta      = truth_muon_eta->at(muon_ind);
        cur_muon.truth_phi      = truth_muon_phi->at(muon_ind);
        cur_muon.truth_charge   = truth_muon_ch->at(muon_ind);
        cur_muon.truth_bar      = truth_muon_barcode->at(muon_ind);

        // -------- truth-to-reco-muon matching --------
        int reco_ind = -1;
        auto it = std::find(real_muon_truth_barcode_list.begin(), real_muon_truth_barcode_list.end(), cur_muon.truth_bar);

        if (it != real_muon_truth_barcode_list.end()) { // found reco muon
            cur_muon.reco_match = true;
            reco_ind = std::distance(real_muon_truth_barcode_list.begin(), it);

            cur_muon.pt = fabs(muon_pt->at(muon_ind))/1000.0;
            cur_muon.eta = muon_eta->at(muon_ind);
            cur_muon.phi = muon_phi->at(muon_ind);
            cur_muon.charge = (muon_pt ->at(muon_ind) > 0)? 1:-1;
            
            cur_muon.dP_overP = muon_deltaP_overP->at(muon_ind);
            cur_muon.z0 = muon_z0->at(muon_ind);
            cur_muon.d0 = muon_d0->at(muon_ind);
            cur_muon.quality = muon_quality->at(muon_ind);
            
            cur_muon.trk_pt      = fabs(muon_trk_pt->at(muon_ind))/1000.0;
            cur_muon.trk_eta     = muon_trk_eta->at(muon_ind);
            cur_muon.trk_phi     = muon_trk_phi->at(muon_ind);

            if (turn_on_track_charge){
                cur_muon.trk_charge = (muon_trk_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge
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

        if(self()->output_single_muon_tree) self()->FillSingleMuonTree();
        
        truth_muon_list.push_back(std::move(cur_muon));
    } // end loop over truth muons

    if(self()->output_single_muon_tree) return;

    // -------- build truth pairs --------

    for (int i = 0; i < truth_muon_list.size(); i++){
        for (int j = 0; j < truth_muon_list.size(); j++){
            mpair().Clear();

            mpair()->Q          = Q;
            mpair()->weight     = static_cast<double>(EventWeights->at(0) * filter_effcy);
            mpair()->crossx     = EventWeights->at(0);

            mpair()->m1 = truth_muon_list.at(i); // default copy constructor
            mpair()->m2 = truth_muon_list.at(j); // default copy constructor

            mpair()->pair_pass_medium = (mpair()->m1.pass_medium && mpair()->m2.pass_medium);
            mpair()->pair_pass_tight  = (mpair()->m1.pass_tight  && mpair()->m2.pass_tight);
        }
    }
    
    self()->FillMuonPairTree();
}
