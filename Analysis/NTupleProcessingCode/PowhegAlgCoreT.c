#include "PowhegAlgCoreT.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"
#include <math.h> 
#include <assert.h>

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::InitParams_PowhegCore(){
    perform_truth |= (!is_fullsim);
    this->output_single_muon_tree &= is_fullsim; // should not output single-muon trees for MC truth

    this->cutLabels = cutLabels_MC;
    this->numCuts = cutLabels_MC.size();

    std::cout << "file_batch: " << file_batch << std::endl;
    std::cout << "is_fullsim? " << is_fullsim << std::endl;
    std::cout << "is_fullsim_overlay? " << is_fullsim_overlay << std::endl;
    std::cout << "perform_truth? " << perform_truth << std::endl;


    mcdir = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/";
    data_subdir = is_fullsim ? "user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim." + mc_mode + ".June2024.1._MYSTREAM/"
                             : mc_mode + "_evgen_truth_full_sample/";

    dt_suffix = is_fullsim? "_fullsim" : "_truth";

    truth_suffix = is_fullsim   ?   (perform_truth ? "_w_truth" : "_no_truth")
                                            :   "";

    if (!is_fullsim && file_batch > 6){
        throw std::runtime_error("For Powheg MC truth, file_batch has to be between 1 and 6!");
    }

    if (is_fullsim && file_batch > 11){
        throw std::runtime_error("For Powheg MC truth, file_batch has to be between 1 and 6!");
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

    fChain() = new TChain("HeavyIonD3PD","HeavyIonD3PD");
    fChain()->SetMakeClass(1);

    for (int ifile = 5 * (file_batch - 1); ifile < 5 * file_batch; ifile++){
        char filename[100];
        if (is_fullsim)        std::sprintf(filename, "%s%suser.yuhang.39654549.MYSTREAM._%06d.root", mcdir.c_str(), data_subdir.c_str(), ifile+1);
        else                   std::sprintf(filename, "%s%smc_truth_%s_%02d.root", mcdir.c_str(), data_subdir.c_str(), mc_mode.c_str(), ifile+1);
        std::cout << filename << std::endl;
        std::ifstream in_file(filename);
        
        if (!in_file.good()){
            std::cout << "Warning: File " << filename << " not found. Skip.\n";
            continue; // skip this file
        }
        fChain()->Add(filename);
    }
   
    cout << "nentries: " << fChain()->GetEntries() << endl;

    fChain()->SetBranchAddress("EventWeights"               , &EventWeights);
    fChain()->SetBranchAddress("Q"                          , &Q);

    fChain()->SetBranchAddress("truth_mupair_pt1"           , &truth_mupair_pt1);
    fChain()->SetBranchAddress("truth_mupair_eta1"          , &truth_mupair_eta1);
    fChain()->SetBranchAddress("truth_mupair_phi1"          , &truth_mupair_phi1);
    fChain()->SetBranchAddress("truth_mupair_ch1"           , &truth_mupair_ch1);
    fChain()->SetBranchAddress("truth_mupair_bar1"          , &truth_mupair_bar1);

    fChain()->SetBranchAddress("truth_mupair_pt2"           , &truth_mupair_pt2);
    fChain()->SetBranchAddress("truth_mupair_eta2"          , &truth_mupair_eta2);
    fChain()->SetBranchAddress("truth_mupair_phi2"          , &truth_mupair_phi2);
    fChain()->SetBranchAddress("truth_mupair_ch2"           , &truth_mupair_ch2);
    fChain()->SetBranchAddress("truth_mupair_bar2"          , &truth_mupair_bar2);

    fChain()->SetBranchAddress("truth_mupair_asym"          , &truth_mupair_asym);
    fChain()->SetBranchAddress("truth_mupair_acop"          , &truth_mupair_acop);
    fChain()->SetBranchAddress("truth_mupair_pt"            , &truth_mupair_pt);
    fChain()->SetBranchAddress("truth_mupair_y"             , &truth_mupair_y);
    fChain()->SetBranchAddress("truth_mupair_m"             , &truth_mupair_m);
  
    fChain()->SetBranchStatus("*"                             ,0);//switch off all branches, then enable just the ones that we need
    fChain()->SetBranchStatus("EventWeights"                  ,1);
    fChain()->SetBranchStatus("Q"                             ,1);

    fChain()->SetBranchStatus("truth_mupair_pt1"              ,1);
    fChain()->SetBranchStatus("truth_mupair_eta1"             ,1);
    fChain()->SetBranchStatus("truth_mupair_phi1"             ,1);
    fChain()->SetBranchStatus("truth_mupair_ch1"              ,1);
    fChain()->SetBranchStatus("truth_mupair_bar1"             ,1);

    fChain()->SetBranchStatus("truth_mupair_pt2"              ,1);
    fChain()->SetBranchStatus("truth_mupair_eta2"             ,1);
    fChain()->SetBranchStatus("truth_mupair_phi2"             ,1);
    fChain()->SetBranchStatus("truth_mupair_ch2"              ,1);
    fChain()->SetBranchStatus("truth_mupair_bar2"             ,1);

    fChain()->SetBranchStatus("truth_mupair_asym"             ,1);
    fChain()->SetBranchStatus("truth_mupair_acop"             ,1);
    fChain()->SetBranchStatus("truth_mupair_pt"               ,1);
    fChain()->SetBranchStatus("truth_mupair_y"                ,1);
    fChain()->SetBranchStatus("truth_mupair_m"                ,1);   

}


template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::OutputTreePathHook(){
    this->output_file_path = mcdir + data_subdir;
    this->output_file_path += "muon_pairs_powheg" + dt_suffix + "_" + mc_mode + "_part" + std::to_string(file_batch) + truth_suffix + ".root";
}

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::InitOutputTreesExtra_PowhegCore(){
    meta_tree = new TTree("meta_tree","meta_tree");
    meta_tree->Branch("nentries_before_cuts", &this->nentries, "nentries_before_cuts/L");    
}

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::OutputHistPathHook(){
    this->output_hist_file_path = mcdir + data_subdir;
    this->output_hist_file_path += "hists_powheg_ntuple_processing_powheg" + dt_suffix + "_" + mc_mode + "_part" + std::to_string(file_batch) + truth_suffix + ".root";
}

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::FillMuonPair_PowhegCore(int pair_ind){
    mpair()->Q          = Q;
    // mpair()->weight     = static_cast<double>(EventWeights->at(0) * filter_effcy / this->nentries);
    mpair()->weight     = static_cast<double>(EventWeights->at(0) * filter_effcy);
    mpair()->crossx     = EventWeights->at(0);

    mpair()->m1.truth_bar     = truth_mupair_bar1->at(pair_ind);
    mpair()->m2.truth_bar     = truth_mupair_bar2->at(pair_ind);
    mpair()->m1.truth_charge  = truth_mupair_ch1->at(pair_ind);//sign of pt stores charge
    mpair()->m2.truth_charge  = truth_mupair_ch2->at(pair_ind);//sign of pt stores charge

    mpair()->m1.truth_pt    = fabs(truth_mupair_pt1->at(pair_ind))/1000.0;//pt of the first muon in the pair
    mpair()->m2.truth_pt    = fabs(truth_mupair_pt2->at(pair_ind))/1000.0;//pt of the second muon in the pair
    mpair()->m1.truth_eta   = truth_mupair_eta1->at(pair_ind);
    mpair()->m2.truth_eta   = truth_mupair_eta2->at(pair_ind);
    mpair()->m1.truth_phi   = truth_mupair_phi1->at(pair_ind);
    mpair()->m2.truth_phi   = truth_mupair_phi2->at(pair_ind);

    mpair()->truth_minv       = truth_mupair_m->at(pair_ind)/1000.;
    mpair()->truth_pair_pt    = truth_mupair_pt->at(pair_ind)/1000.;
    mpair()->truth_pair_y     = truth_mupair_y->at(pair_ind);
    mpair()->truth_asym       = abs(truth_mupair_asym->at(pair_ind));
    mpair()->truth_acop       = abs(truth_mupair_acop->at(pair_ind));
    assert (mpair()->truth_asym >= 0 && mpair()->truth_acop >= 0);
}

template <class PairT, class MuonT, class Derived, class... Extras>
bool PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::PassCuts_PowhegCore(){
    // truth is the base (reco is with detector effect --> relative to truth): apply cuts w.r.t. truth kinematics
    if (fabs(mpair()->m1.truth_eta) > 2.4 || fabs(mpair()->m2.truth_eta) > 2.4) return false;
    h_cutAcceptance()[mpair()->m1.truth_charge != mpair()->m2.truth_charge]->Fill(double(pass_muon_eta) + 0.5, mpair()->weight); // if same sign: fill the h_cutAcceptance()[0] histogram; if opposite sign, fill the h_cutAcceptance()[1] histogram

    if (mpair()->m1.truth_pt < 4 || mpair()->m2.truth_pt < 4) return false;
    h_cutAcceptance()[mpair()->m1.truth_charge != mpair()->m2.truth_charge]->Fill(double(pass_muon_pt) + 0.5, mpair()->weight);

    return true;
 }

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras...>::ProcessDataHook(){
    this->nentries = fChain()->GetEntries();//number of events
    meta_tree->Fill();
    
    // for (Long64_t jentry=0; jentry<1000; jentry++) {//loop over the events
    for (Long64_t jentry=0; jentry<this->nentries; jentry++) {//loop over the events
        // cout << jentry << endl;
        if(jentry%10000==0) cout<<"Processing "<<jentry<<" event out of "<<this->nentries<<" events"<<std::endl;

        int num_bytes = fChain()->GetEntry(jentry);//read in an event
        if(num_bytes==0){
          std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
          throw std::exception();
        }
    
        this->muon_pair_list_cur_event_pre_resonance_cut.clear();
        this->resonance_tagged_muon_index_list.clear(); // MUST CLEAR for each event!!

        std::vector<int> muon_index_list = {};
        std::vector<int>::iterator it;

        int NPairs = muon_pair_muon1_pt->size();//number of muon pairs in the event

        for(int pair_ind = 0; pair_ind < NPairs; pair_ind++){//first loop over all muon-pairs in the event

            mpair().Clear();

            this->FillMuonPair(pair_ind);
            
            mpair()->m1.ev_num = jentry;
            mpair()->m2.ev_num = jentry;

            if (perform_truth){ // order: truth > reco
                h_cutAcceptance()[mpair()->m1.truth_charge != mpair()->m2.truth_charge]->Fill(double(nocut) + 0.5, mpair()->weight);
            } else{
                h_cutAcceptance()[mpair()->m1.charge != mpair()->m2.charge]->Fill(double(nocut) + 0.5, mpair()->weight);
            }

            // ------------------------------------------------------------

            // Apply cuts

            // Trigger match for muon pair
            // if(dimuon_b_HLT_2mu4->at(pair_ind)==false) continue;

            if (abs(mpair()->weight) > crossx_cut * filter_effcy) continue;
            if (!PassCuts(mpair()))continue;
        
            //------------------------------------------------------------

            mpair()->Update(); // only afterwards can we use mpair()->same_sign

            // resonance tag
            ResonanceTagging(mpair());
            
            this->muon_pair_list_cur_event_pre_resonance_cut.push_back(std::move(mpair()));
            
        } // finish first loop over all muon pairs

        for(int pair_ind = 0; pair_ind < this->muon_pair_list_cur_event_pre_resonance_cut.size(); pair_ind++){//second loop over all muon-pairs in the event
            // discard the pair if either muon is resonance-tagged

            mpair() = std::move(this->muon_pair_list_cur_event_pre_resonance_cut.at(pair_ind));

            if (!mpair()){
                std::cerr << "mpair() at second muon-pair loop NOT found! Return without resonance-tag checking or pair analysis." << std::endl;
                continue;
            }

            std::vector<int>::iterator itres_m1;
            std::vector<int>::iterator itres_m2;

            itres_m1 = std::find(this->resonance_tagged_muon_index_list.begin(),this->resonance_tagged_muon_index_list.end(),mpair()->m1.ind);
            if(itres_m1 != this->resonance_tagged_muon_index_list.end())  continue;

            itres_m2 = std::find(this->resonance_tagged_muon_index_list.begin(),this->resonance_tagged_muon_index_list.end(),mpair()->m2.ind);
            if(itres_m2 != this->resonance_tagged_muon_index_list.end())  continue;

            h_cutAcceptance()[mpair()->m1.charge != mpair()->m2.charge]->Fill(double(pass_resonance) + 0.5, mpair()->weight);
            
            //------------------------------------------------------------

            if (is_fullsim){
                if(this->output_single_muon_tree){
                  it = std::find(muon_index_list.begin(),muon_index_list.end(),mpair()->m1.ind);
                  if(it == muon_index_list.end()){ //muon1 index NOT found
                        muon_index_list.push_back(mpair()->m1.ind);
                        this->muon_raw_ptr = &(mpair()->m1);
                        this->FillSingleMuonTree();
                  }
                  it = std::find(muon_index_list.begin(),muon_index_list.end(),mpair()->m2.ind);
                  if(it == muon_index_list.end()){ //muon1 index NOT found
                        muon_index_list.push_back(mpair()->m2.ind);
                        this->muon_raw_ptr = &(mpair()->m2);
                        this->FillSingleMuonTree();
                  }
                }
                else{ // fill muon pair trees
                    this->FillMuonPairTree();
                }
            }

            PerformTruthPairAnalysisHook();

        } // finish second loop over muon pairs

        this->muon_pair_list_cur_event_pre_resonance_cut.clear();
    }//loop over events
}

