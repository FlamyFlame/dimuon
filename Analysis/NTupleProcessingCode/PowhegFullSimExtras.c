#include "PowhegFullSimExtras.h"
#include "../Utilities/tchain_helpers.h"
#include "../Utilities/AllVertexIPSelection.h"

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

    // ALL-VERTEX impact-parameter selection (pp-conditions fullsim only). enable_and_bind
    // enables the branch BEFORE setting its address, which is the order that matters on a
    // SetMakeClass(1) chain -- binding a still-DISABLED STL-collection branch after the first
    // tree is loaded leaves the pointer null. See Utilities/AllVertexIPSelection.h
    // (BRANCH BINDING).
    if (UseAllVertexIP()) {
        if (!self().fChainRef()->GetBranch("vtx_z") ||
            !self().fChainRef()->GetBranch("vtx_ntrk"))
            throw std::runtime_error("PowhegFullSimExtras: pp-conditions fullsim NTUP has no "
                                     "vtx_z/vtx_ntrk branches, so the all-vertex "
                                     "impact-parameter selection that MIRRORS the pp data "
                                     "selection cannot be applied. Re-skim with m_store_Vtx.");
        enable_and_bind(self().fChainRef(), "vtx_z"   , &vtx_z);
        enable_and_bind(self().fChainRef(), "vtx_ntrk", &vtx_ntrk);
    }
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
    // 2026-09-08: fabs restored. Without it this kept muons with eta < -2.4, i.e. it was NOT the
    // mirror of PassCuts_DataCore that the comment above claims (Pythia's twin always had the
    // fabs). Listed as an open defect in muon_gap_cuts_acceptance.md F17; closed here because it
    // sits inside the block being edited. Only reachable on the POWHEG fullsim RECO path, which
    // has no live consumer -- see pp24_all_vertex_pairs.md ("Why POWHEG is code-changed but NOT
    // rerun"); the trees on disk still carry the un-fabs'd selection.
    if (fabs(muon.eta) > 2.4) return false;
    // Reco threshold 4.0 -> 4.5 GeV (user decision 2026-09-08,
    // mu_pt45_gap125_pairpt9_adoption.md §3(a)); mirrors DimuonDataAlgCoreT::PassCuts_DataCore.
    if (muon.pt < 4.5) return false;

    // HF muon cut    
    if (muon.dP_overP > self().pmsRef().deltaP_overP_thrsh ) return false;
    
    //cut on d0 & z0 sin(theta) against fake muons
    if (UseAllVertexIP()) {
        // ANY-vertex form for the per-muon flag; the stricter SAME-vertex requirement is added
        // at pair level in PassPairAllVertexIP (see PythiaFullSimExtras.h).
        if (!AllVertexIP::PassAnyVertex(muon.d0, muon.z0, muon.eta, vtx_z, vtx_ntrk,
                                        self().pmsRef().d0cut, self().pmsRef().z0cut))
            return false;
    } else {
        double z0sinTheta = fabs(muon.z0 * sin(2.0*atan(exp(-muon.eta))));
        bool pass_d0_z0_cuts = (fabs(muon.d0) < self().pmsRef().d0cut && z0sinTheta < self().pmsRef().z0cut);
        if (!pass_d0_z0_cuts) return false;
    }

    // if track charge saved, require muon + track charge to agree
    if (turn_on_track_charge){
        if (muon.trk_charge != muon.charge) return false;
    }
    
    return true;
}

template <class PairT, class MuonT, class Derived>
int PowhegFullSimExtras<PairT, MuonT, Derived>::PairAllVertexIndex(const muon_t& m1, const muon_t& m2,
                                                                    bool* pass_primary_out){
    // Pure lookup, NO counters -- called before the WP flags are known.
    return AllVertexIP::BestCommonVertex(m1.d0, m1.z0, m1.eta,
                                         m2.d0, m2.z0, m2.eta,
                                         vtx_z, vtx_ntrk,
                                         self().pmsRef().d0cut, self().pmsRef().z0cut,
                                         pass_primary_out);
}


template <class PairT, class MuonT, class Derived>
void PowhegFullSimExtras<PairT, MuonT, Derived>::FinalizeExtra(){
    if (!UseAllVertexIP() || n_allvtx_pairs_pass == 0) return;
    auto pct = [](long long a, long long b){ return b > 0 ? 100.0 * a / b : 0.0; };
    std::cout << "All-vertex IP report (POWHEG pp fullsim; both muons reco-matched + Tight WP; "
                 "no trigger, no resonance veto): "
              << n_allvtx_pairs_pass << " pairs\n"
              << "  from a SECONDARY vertex (best vertex != 0): " << n_allvtx_pairs_secondary
              << " (" << pct(n_allvtx_pairs_secondary, n_allvtx_pairs_pass) << "%)\n"
              << "  ADDED by the all-vertex rule (fail primary): " << n_allvtx_pairs_added
              << " (" << pct(n_allvtx_pairs_added, n_allvtx_pairs_pass) << "%)" << std::endl;
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
            // Mirror the pair-level truth cuts from PassCuts_PowhegCore():
            // only store muons that are in the fiducial acceptance used for
            // the muon-pair analysis (truth_pt > 4.5 GeV, |truth_eta| < 2.4).
            // Threshold 4.0 -> 4.5 GeV (user decision 2026-09-08,
            // mu_pt45_gap125_pairpt9_adoption.md §3(a)).
            if (cur_muon.truth_pt > 4.5 && fabs(cur_muon.truth_eta) < 2.4){
                self().muon_raw_ptr = &cur_muon;
                self().FillSingleMuonTree();
            }
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
                // Same-vertex requirement (pp-conditions fullsim only), mirroring the pp data.
                bool ip_pass_primary = false;
                const int ip_vtx = UseAllVertexIP()
                    ? PairAllVertexIndex(self().mpairRef()->m1, self().mpairRef()->m2, &ip_pass_primary)
                    : -1;
                // No `disable_ip_cut` term here, unlike the Pythia mirror: PowhegFullSimExtras
                // has no such member. The asymmetry is deliberate, not an oversight.
                const bool ip_pair_ok = !UseAllVertexIP() || ip_vtx >= 0;
                self().mpairRef()->pair_pass_medium = (self().mpairRef()->m1.pass_medium && self().mpairRef()->m2.pass_medium && ip_pair_ok);
                self().mpairRef()->pair_pass_tight  = (self().mpairRef()->m1.pass_tight  && self().mpairRef()->m2.pass_tight  && ip_pair_ok);

                // Counted only over Tight-WP pairs, so the reported fraction has a stated
                // population (see FinalizeExtra).
                if (UseAllVertexIP() && ip_vtx >= 0 && self().mpairRef()->pair_pass_tight) {
                    ++n_allvtx_pairs_pass;
                    if (ip_vtx != 0)      ++n_allvtx_pairs_secondary;
                    if (!ip_pass_primary) ++n_allvtx_pairs_added;
                }
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
