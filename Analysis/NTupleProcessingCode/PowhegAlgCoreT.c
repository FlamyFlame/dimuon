#include "PowhegAlgCoreT.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"
#include <math.h> 
#include <assert.h>
#include "../Utilities/tchain_helpers.h"


template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::InitParams_PowhegCore(){
    this->cutLabels = cutLabels_MC;
    this->numCuts = cutLabels_MC.size();
    this->isMC = true; // is MC: regardless of truth only or fullsim

    std::cout << "file_batch: " << file_batch << std::endl;
    std::cout << "is_fullsim? " << is_fullsim << std::endl;
    std::cout << "is_fullsim_overlay? " << is_fullsim_overlay << std::endl;
    std::cout << "perform_truth? " << perform_truth << std::endl;

    if (is_fullsim && (run_year < 15 || run_year > 26)){
        throw std::runtime_error("For Powheg MC fullsim/fullsim overlay, run_year has to be between 15 and 26.");
    }

    powheg_dir = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/";
    
    std::string run_year_str = std::to_string(run_year);

    if (is_fullsim_overlay){
        data_subdir = "user.yuhang.TrigRates.dimuon.PowhegPythia.fullsimOverlay.PbPb" + run_year_str + "." + mc_mode + ".Feb2026.v1._MYSTREAM/";
    } else if (is_fullsim){
        data_subdir = "user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp" + run_year_str + "." + mc_mode + ".Feb2026.v1._MYSTREAM/";
    } else{
        data_subdir = mc_mode + "_evgen_truth_full_sample/";
    }

    dt_suffix = "";
    if (is_fullsim_overlay){
        dt_suffix = (perform_truth ? "_fullsim_overlay_PbPb" + run_year_str + "_w_truth" : "_fullsim_overlay_PbPb" + run_year_str + "_no_truth");
    } else if (is_fullsim){
        dt_suffix = (perform_truth ? "_fullsim_pp" + run_year_str + "" : "_fullsim_pp" + run_year_str + "_no_truth");
    } else{
        dt_suffix = "_truth";
    }

    if (!is_fullsim && file_batch > 6){ // truth
        throw std::runtime_error("For Powheg MC truth, file_batch has to be between 1 and 6!");
    }

    if (is_fullsim && (!is_fullsim_overlay) && file_batch > 11) { // fullsim, no overlay
        throw std::runtime_error("For Powheg MC fullsim, file_batch has to be between 1 and 11!");
    }

    if (is_fullsim_overlay && file_batch > 11){ // fullsim overlay
        throw std::runtime_error("For Powheg MC fullsim overlay, file_batch has to be between 1 and 11!");
    }

    if (mc_mode == "bb") filter_effcy = filter_effcy_bb;
    else if (mc_mode == "cc") filter_effcy = filter_effcy_cc;
    else{
        filter_effcy = -1000;
        std::cout << "MC mode has to be 'bb' or 'cc.' Else program will terminate." << std::endl;
        throw std::exception();
    }
    std::cout << "The MC mode is " << mc_mode << ". Filter efficiency is " << filter_effcy << std::endl;
}

//initialize the TChain
template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::InitInput_PowhegCore(){

    fChainRef() = new TChain("HeavyIonD3PD","HeavyIonD3PD");

    fChainRef()->SetMakeClass(1);

    std::string task_id = (mc_mode == "bb")? "48591366" : "48591379";

    for (int ifile = 5 * (file_batch - 1); ifile < 5 * file_batch; ifile++){
        std::string filename =
            is_fullsim
            ? (powheg_dir + data_subdir + "user.yuhang." + task_id + ".MYSTREAM._" + Form("%06d", ifile+1) + ".root")
            : (powheg_dir + data_subdir + "mc_truth_" + mc_mode + "_" + Form("%02d", ifile+1) + ".root");

        std::cout << filename << std::endl;
        std::ifstream in_file(filename.c_str());
        
        if (!in_file.good()){
            std::cout << "Warning: File " << filename << " not found. Skip.\n";
            continue; // skip this file
        }
        fChainRef()->Add(filename.c_str());
    }
    
    fChainRef()->LoadTree(0);
    cout << "nentries: " << fChainRef()->GetEntries() << endl;

    fChainRef()->SetBranchStatus("*"                             ,0);//switch off all branches, then enable just the ones that we need

    enable_and_bind(fChainRef(), "EventWeights"               , &EventWeights);

    enable_and_bind(fChainRef(), "truth_muon_barcode"         , &truth_muon_barcode);

    enable_and_bind(fChainRef(), "truth_mupair_pt1"           , &truth_mupair_pt1);
    enable_and_bind(fChainRef(), "truth_mupair_eta1"          , &truth_mupair_eta1);
    enable_and_bind(fChainRef(), "truth_mupair_phi1"          , &truth_mupair_phi1);
    enable_and_bind(fChainRef(), "truth_mupair_ch1"           , &truth_mupair_ch1);
    enable_and_bind(fChainRef(), "truth_mupair_bar1"          , &truth_mupair_bar1);

    enable_and_bind(fChainRef(), "truth_mupair_pt2"           , &truth_mupair_pt2);
    enable_and_bind(fChainRef(), "truth_mupair_eta2"          , &truth_mupair_eta2);
    enable_and_bind(fChainRef(), "truth_mupair_phi2"          , &truth_mupair_phi2);
    enable_and_bind(fChainRef(), "truth_mupair_ch2"           , &truth_mupair_ch2);
    enable_and_bind(fChainRef(), "truth_mupair_bar2"          , &truth_mupair_bar2);

    enable_and_bind(fChainRef(), "truth_mupair_asym"          , &truth_mupair_asym);
    enable_and_bind(fChainRef(), "truth_mupair_acop"          , &truth_mupair_acop);
    enable_and_bind(fChainRef(), "truth_mupair_pt"            , &truth_mupair_pt);
    enable_and_bind(fChainRef(), "truth_mupair_y"             , &truth_mupair_y);
    enable_and_bind(fChainRef(), "truth_mupair_m"             , &truth_mupair_m);
}

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::CheckBranchPtrs_PowhegCore(){
    if (this->debug_mode) std::cout << "Calling CheckBranchPtrs_PowhegCore" << std::endl;
    
    auto require = [&](auto* p, const char* name){
        if(!p) throw std::runtime_error(std::string("Null branch pointer: ") + name);
    };

    require(EventWeights, "EventWeights");
    require(truth_muon_barcode, "truth_muon_barcode");
    require(truth_mupair_pt1, "truth_mupair_pt1");
    require(truth_mupair_eta1, "truth_mupair_eta1");
    require(truth_mupair_phi1, "truth_mupair_phi1");
    require(truth_mupair_ch1, "truth_mupair_ch1");
    require(truth_mupair_bar1, "truth_mupair_bar1");
    require(truth_mupair_pt2, "truth_mupair_pt2");
    require(truth_mupair_eta2, "truth_mupair_eta2");
    require(truth_mupair_phi2, "truth_mupair_phi2");
    require(truth_mupair_ch2, "truth_mupair_ch2");
    require(truth_mupair_bar2, "truth_mupair_bar2");
    require(truth_mupair_asym, "truth_mupair_asym");
    require(truth_mupair_acop, "truth_mupair_acop");
    require(truth_mupair_pt, "truth_mupair_pt");
    require(truth_mupair_y, "truth_mupair_y");
    require(truth_mupair_m, "truth_mupair_m");

    if (this->debug_mode) std::cout << "truth_mupair_pt1 size " << truth_mupair_pt1->size() << ", " << truth_mupair_pt1->at(0) << std::endl;
    if (this->debug_mode) std::cout << "truth_mupair_bar2 size " << truth_mupair_bar2->size() << ", " << truth_mupair_bar2->at(0) << std::endl;

}


template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::OutputTreePathHook(){
    this->output_file_path = powheg_dir + data_subdir;
    std::string file_name_base = this->output_single_muon_tree ? "single_muon_trees_powheg" : "muon_pairs_powheg";

    this->output_file_path += file_name_base + "_" + mc_mode + dt_suffix + "_part" + std::to_string(file_batch) + ".root";
}

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::InitOutputTreesExtra_PowhegCore(){
    meta_tree = new TTree("meta_tree","meta_tree");
    meta_tree->Branch("nentries_before_cuts", &this->nentries, "nentries_before_cuts/L");    
}

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::OutputHistPathHook(){
    this->output_hist_file_path = powheg_dir + data_subdir;
    this->output_hist_file_path += "hists_powheg_ntuple_processing_powheg_" + mc_mode + dt_suffix + "_part" + std::to_string(file_batch) + ".root";
}

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::FillMuonPair_PowhegCore(int pair_ind){
    float event_crossx = (EventWeights->size() == 0)? 1. : EventWeights->at(0);
    float event_weight = (EventWeights->size() == 0)? 1. : event_crossx * filter_effcy;

    mpairRef()->weight     = event_weight;
    mpairRef()->crossx     = event_crossx;

    mpairRef()->m1.truth_bar     = truth_mupair_bar1->at(pair_ind);
    mpairRef()->m2.truth_bar     = truth_mupair_bar2->at(pair_ind);
    mpairRef()->m1.truth_charge  = truth_mupair_ch1->at(pair_ind);//sign of pt stores charge
    mpairRef()->m2.truth_charge  = truth_mupair_ch2->at(pair_ind);//sign of pt stores charge

    mpairRef()->m1.truth_pt    = fabs(truth_mupair_pt1->at(pair_ind))/1000.0;//pt of the first muon in the pair
    mpairRef()->m2.truth_pt    = fabs(truth_mupair_pt2->at(pair_ind))/1000.0;//pt of the second muon in the pair
    mpairRef()->m1.truth_eta   = truth_mupair_eta1->at(pair_ind);
    mpairRef()->m2.truth_eta   = truth_mupair_eta2->at(pair_ind);
    mpairRef()->m1.truth_phi   = truth_mupair_phi1->at(pair_ind);
    mpairRef()->m2.truth_phi   = truth_mupair_phi2->at(pair_ind);

    mpairRef()->truth_minv       = truth_mupair_m->at(pair_ind)/1000.;
    mpairRef()->truth_pair_pt    = truth_mupair_pt->at(pair_ind)/1000.;
    mpairRef()->truth_pair_y     = truth_mupair_y->at(pair_ind);
    mpairRef()->truth_asym       = abs(truth_mupair_asym->at(pair_ind));
    mpairRef()->truth_acop       = abs(truth_mupair_acop->at(pair_ind));
    assert (mpairRef()->truth_asym >= 0 && mpairRef()->truth_acop >= 0);

    auto it1 = std::find(truth_muon_barcode->begin(), truth_muon_barcode->end(), mpairRef()->m1.truth_bar);
    auto it2 = std::find(truth_muon_barcode->begin(), truth_muon_barcode->end(), mpairRef()->m2.truth_bar);

    if (it1 == truth_muon_barcode->end() || it2 == truth_muon_barcode->end()) { // either muon barcode not matched
        throw std::runtime_error("Truth muon in current pair NOT matched to truth muon list!");
    }
    
    mpairRef()->m1.ind = std::distance(truth_muon_barcode->begin(), it1);
    mpairRef()->m2.ind = std::distance(truth_muon_barcode->begin(), it2);

}

template <class PairT, class MuonT, class Derived, class... Extras>
bool PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::PassCuts_PowhegCore(){
    // truth is the base (reco is with detector effect --> relative to truth): apply cuts w.r.t. truth kinematics
    if (fabs(mpairRef()->m1.truth_eta) > 2.4 || fabs(mpairRef()->m2.truth_eta) > 2.4) return false;
    h_cutAcceptanceRef()[mpairRef()->m1.truth_charge != mpairRef()->m2.truth_charge]->Fill(double(pass_muon_eta) + 0.5, mpairRef()->weight); // if same sign: fill the h_cutAcceptanceRef()[0] histogram; if opposite sign, fill the h_cutAcceptanceRef()[1] histogram

    if (mpairRef()->m1.truth_pt < 4 || mpairRef()->m2.truth_pt < 4) return false;
    h_cutAcceptanceRef()[mpairRef()->m1.truth_charge != mpairRef()->m2.truth_charge]->Fill(double(pass_muon_pt) + 0.5, mpairRef()->weight);

    return true;
 }

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::ProcessDataHook(){
    this->nentries = fChainRef()->GetEntries();//number of events

    const Long64_t nentries_to_process = (this->nevents_max <= 0)
        ? this->nentries
        : std::min<Long64_t>(this->nentries, this->nevents_max);

    Long64_t nentries_actual = nentries_to_process;
    
    // for (Long64_t jentry=0; jentry<10000; jentry++) {//loop over the events
    for (Long64_t jentry=0; jentry<nentries_to_process; jentry++) {//loop over the events
        if (this->debug_mode)       std::cout << "Event#: " << jentry << std::endl;
        else if(jentry%10000==0)    std::cout << "Processing "<<jentry<<" event out of "<<nentries_to_process<<" events"<<std::endl;

        int num_bytes = fChainRef()->GetEntry(jentry);//read in an event
        if(num_bytes<=0){
          std::cerr<<"Error:: Read in event has size of zero bytes,  skipping"<<std::endl;
          nentries_actual--;
          continue;
        }

        CheckBranchPtrs();

        if (!is_fullsim)    ProcessEventTruthOnly(jentry);
        else                ProcessEventFullsimHook(jentry);
    }

    std::cout << "#Entries in file: " << this->nentries << ", #Entries processed: " << nentries_to_process << ", #Entries with non-zero size: " << nentries_actual << std::endl;
    this->nentries = nentries_actual;
    meta_tree->Fill();

}

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::ProcessEventTruthOnly(int ev_num){

    this->muon_pair_list_cur_event_pre_resonance_cut.clear();
    this->resonance_tagged_muon_index_list.clear(); // MUST CLEAR for each event!!

    int NPairs = truth_mupair_pt1->size();//number of muon pairs in the event

    for(int pair_ind = 0; pair_ind < NPairs; pair_ind++){//first loop over all muon-pairs in the event

        mpairRef() = std::make_shared<pair_t>();

        this->FillMuonPair(pair_ind);
        
        mpairRef()->m1.ev_num = ev_num;
        mpairRef()->m2.ev_num = ev_num;

        h_cutAcceptanceRef()[mpairRef()->m1.truth_charge != mpairRef()->m2.truth_charge]->Fill(double(nocut) + 0.5, mpairRef()->weight);

        // ------------------------------------------------------------

        // Apply cuts

        if (abs(mpairRef()->weight) > crossx_cut * filter_effcy) continue;
        if (!this->PassCuts())continue;
        
        //------------------------------------------------------------

        mpairRef()->Update(); // only afterwards can we use mpairRef()->same_sign

        // resonance tag
        this->ResonanceTagging();
        
        this->muon_pair_list_cur_event_pre_resonance_cut.push_back(std::move(mpairRef()));
        
    } // finish first loop over all muon pairs

    for(int pair_ind = 0; pair_ind < this->muon_pair_list_cur_event_pre_resonance_cut.size(); pair_ind++){//second loop over all muon-pairs in the event
        // discard the pair if either muon is resonance-tagged

        mpairRef() = std::move(this->muon_pair_list_cur_event_pre_resonance_cut.at(pair_ind));

        if (!mpairRef()){
            std::cerr << "mpairRef() at second muon-pair loop NOT found! Skipping pair in output tree filling." << std::endl;
            continue;
        }

        std::vector<int>::iterator itres_m1;
        std::vector<int>::iterator itres_m2;

        itres_m1 = std::find(this->resonance_tagged_muon_index_list.begin(),this->resonance_tagged_muon_index_list.end(),mpairRef()->m1.ind);
        if(itres_m1 != this->resonance_tagged_muon_index_list.end())  continue;

        itres_m2 = std::find(this->resonance_tagged_muon_index_list.begin(),this->resonance_tagged_muon_index_list.end(),mpairRef()->m2.ind);
        if(itres_m2 != this->resonance_tagged_muon_index_list.end())  continue;

        h_cutAcceptanceRef()[mpairRef()->m1.truth_charge != mpairRef()->m2.truth_charge]->Fill(double(pass_resonance) + 0.5, mpairRef()->weight);
        
        PerformTruthPairAnalysisHook();
        this->FillMuonPairTree();
        
    } // finish second loop over muon pairs

    mpairRef().reset(); // make sure the reference to the last muon pair gets reset
}

