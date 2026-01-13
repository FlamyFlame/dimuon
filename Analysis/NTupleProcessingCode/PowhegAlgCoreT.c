#include "PowhegAlgCoreT.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"
#include <math.h> 
#include <assert.h>

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras>::InitParams_PowhegCore(){
    perform_truth |= (!is_fullsim);
    output_single_muon_tree &= is_fullsim; // should not output single-muon trees for MC truth

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
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras>::InitInput_PowhegCore(){

    fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
    fChain->SetMakeClass(1);

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
        fChain->Add(filename);
    }
   
    cout << "nentries: " << fChain->GetEntries() << endl;

    fChain->SetBranchAddress("EventWeights"               , &EventWeights);
    fChain->SetBranchAddress("Q"               , &Q);

    fChain->SetBranchStatus("*"                               ,0);//switch off all branches, then enable just the ones that we need
    fChain->SetBranchStatus("EventWeights"             ,1);
    fChain->SetBranchStatus("Q"             ,1);

}


template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras>::OutputTreePath(){
    output_file_path = mcdir + data_subdir;
    output_file_path += "muon_pairs_powheg" + dt_suffix + "_" + mc_mode + "_part" + std::to_string(file_batch) + truth_suffix + ".root";
}

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras>::InitOutputTreesExtra_PowhegCore(){
    meta_tree = new TTree("meta_tree","meta_tree");
    meta_tree->Branch("nentries_before_cuts", &nentries, "nentries_before_cuts/L");    
}

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras>::OutputHistPath(){
    output_hist_file_path = mcdir + data_subdir;
    output_hist_file_path += "hists_powheg_ntuple_processing_powheg" + dt_suffix + "_" + mc_mode + "_part" + std::to_string(file_batch) + truth_suffix + ".root";
}

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras>::FillMuonPair_PowhegCore(int pair_ind){
  mpair->Q          = Q;
  // mpair->weight     = static_cast<double>(EventWeights->at(0) * filter_effcy / nentries);
  mpair->weight     = static_cast<double>(EventWeights->at(0) * filter_effcy);
  mpair->crossx     = EventWeights->at(0);
}

template <class PairT, class MuonT, class Derived, class... Extras>
bool PowhegAlgCoreT<PairT, MuonT, Derived, Extras>::PassCuts_PowhegCore(){
  //Apply ALL CUTS but for resonances

  // NO quality cuts on the MC truth muons
  // if((mpair->m1.quality&mpair->m2.quality&1  )==0) return false;//compair->m2ined muon
  // if((mpair->m1.quality&mpair->m2.quality&8  )==0) return false;//Medium muon
  // if((mpair->m1.quality&mpair->m2.quality&32 )==0) return false;//IDCuts
  //if((mpair->m1.quality&mpair->m2.quality&256)==0) return false;//MuonCuts

  if (fabs(mpair->m1.eta) > 2.4 || fabs(mpair->m2.eta) > 2.4) return false;
  
  // passed muon eta cut
  // if (mpair->same_sign) h_cutAcceptance[0]->Fill(pass_muon_eta + 0.5);
  // else h_cutAcceptance[1]->Fill(pass_muon_eta + 0.5);
  h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(double(pass_muon_eta) + 0.5, mpair->weight); // if same sign: fill the h_cutAcceptance[0] histogram; if opposite sign, fill the h_cutAcceptance[1] histogram

  if (mpair->m1.pt < 4 || mpair->m2.pt < 4) return false;
  
  // passed muon pt cut
  // if (mpair->same_sign) h_cutAcceptance[0]->Fill(pass_muon_pt + 0.5);
  // else h_cutAcceptance[1]->Fill(pass_muon_pt + 0.5);
  h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(double(pass_muon_pt) + 0.5, mpair->weight);

  // if( fabs(mpair->m1.dP_overP) > pms.deltaP_overP_thrsh || fabs(mpair->m2.dP_overP) > pms.deltaP_overP_thrsh ) pass = false;
 }

template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras>::ProcessData(){
    nentries = fChain->GetEntries();//number of events
    meta_tree->Fill();
    
    // for (Long64_t jentry=0; jentry<1000; jentry++) {//loop over the events
    for (Long64_t jentry=0; jentry<nentries; jentry++) {//loop over the events
        // cout << jentry << endl;
        if(jentry%10000==0) cout<<"Processing "<<jentry<<" event out of "<<nentries<<" events"<<std::endl;

        int num_bytes = fChain->GetEntry(jentry);//read in an event
        if(num_bytes==0){
          std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
          throw std::exception();
        }
    
        muon_pair_list_cur_event_pre_resonance_cut.clear();
        resonance_tagged_muon_index_list.clear(); // MUST CLEAR for each event!!

        std::vector<int> muon_index_list = {};
        std::vector<int>::iterator it;

        int NPairs = muon_pair_muon1_pt->size();//number of muon pairs in the event

        for(int pair_ind = 0; pair_ind < NPairs; pair_ind++){//first loop over all muon-pairs in the event

            mpair.Clear();

            FillMuonPair(pair_ind);
            
            mpair->m1.ev_num = jentry;
            mpair->m2.ev_num = jentry;

            if (is_fullsim){
                h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(double(nocut) + 0.5, mpair->weight);
            } else{
                h_cutAcceptance[mpair->m1.truth_charge != mpair->m2.truth_charge]->Fill(double(nocut) + 0.5, mpair->weight);
            }

            // ------------------------------------------------------------

            // Apply cuts

            // Trigger match for muon pair
            // if(dimuon_b_HLT_2mu4->at(pair_ind)==false) continue;

            if (abs(mpair->weight) > crossx_cut * filter_effcy) continue;
            if (!PassCuts(mpair))continue;
        
            //------------------------------------------------------------

            mpair->Update(); // only afterwards can we use mpair->same_sign

            // resonance tag
            ResonanceTagging(mpair);

            // photo-production cut - DO NOT APPLY FOR MC
            // if (IsPhotoProduction()) continue;
            
            muon_pair_list_cur_event_pre_resonance_cut.push_back(std::move(mpair));
            
        } // finish first loop over all muon pairs

        for(int pair_ind = 0; pair_ind < muon_pair_list_cur_event_pre_resonance_cut.size(); pair_ind++){//second loop over all muon-pairs in the event
            // discard the pair if either muon is resonance-tagged

            mpair = std::move(muon_pair_list_cur_event_pre_resonance_cut.at(pair_ind));

            if (!mpair){
                std::cerr << "mpair at second muon-pair loop NOT found! Return without resonance-tag checking or pair analysis." << std::endl;
                continue;
            }

            std::vector<int>::iterator itres_m1;
            std::vector<int>::iterator itres_m2;

            itres_m1 = std::find(resonance_tagged_muon_index_list.begin(),resonance_tagged_muon_index_list.end(),mpair->m1.ind);
            if(itres_m1 != resonance_tagged_muon_index_list.end())  continue;

            itres_m2 = std::find(resonance_tagged_muon_index_list.begin(),resonance_tagged_muon_index_list.end(),mpair->m2.ind);
            if(itres_m2 != resonance_tagged_muon_index_list.end())  continue;

            h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(double(pass_resonance) + 0.5, mpair->weight);
            
            //------------------------------------------------------------

            if (is_fullsim){
                if(output_single_muon_tree){
                  it = std::find(muon_index_list.begin(),muon_index_list.end(),mpair->m1.ind);
                  if(it == muon_index_list.end()){ //muon1 index NOT found
                        muon_index_list.push_back(mpair->m1.ind);
                        muon_raw_ptr = &(mpair->m1);
                        FillSingleMuonTree();
                  }
                  it = std::find(muon_index_list.begin(),muon_index_list.end(),mpair->m2.ind);
                  if(it == muon_index_list.end()){ //muon1 index NOT found
                        muon_index_list.push_back(mpair->m2.ind);
                        muon_raw_ptr = &(mpair->m2);
                        FillSingleMuonTree();
                  }
                }
                else{ // fill muon pair trees
                    FillMuonPairTree();
                }
            }

            self().PerformTruthPairAnalysisHook();

        } // finish second loop over muon pairs

        muon_pair_list_cur_event_pre_resonance_cut.clear();
    }//loop over events
}



template <class PairT, class MuonT, class Derived, class... Extras>
void PowhegAlgCoreT<PairT, MuonT, Derived, Extras>::HistAdjust(){
  DimuonAnalysisBaseClass::HistAdjust();
  self().HistAdjustHook();
}
