#ifndef __PowhegNTupleFirstPass_C__
#define __PowhegNTupleFirstPass_C__

#include "PowhegNTupleFirstPass.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"
#include <math.h> 
#include <assert.h>

void PowhegNTupleFirstPass::Run(){
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    Initialize();
    ProcessData();
    Finalize();

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    std::cout << "CPU time used (in seconds) is " << cpu_time_used << std::endl;

}


void PowhegNTupleFirstPass::Initialize(){
    output_single_muon_tree &= is_fullsim; // should not output single-muon trees for MC truth

    mcdir = "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/";
    data_subdir = is_fullsim ? "user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim." + mc_mode + ".June2024.1._MYSTREAM/"
                             : mc_mode + "_evgen_truth_full_sample/";

    dt_suffix = is_fullsim? "_fullsim" : "_truth";

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

    InitInput();
    InitTempVariables();
    InitOutput();
}

//initialize the TChain
void PowhegNTupleFirstPass::InitInput(){

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

    // MC muons have no quality, d0 or z0 recorded
    fChain->SetBranchAddress("truth_id"                   , &truth_id);
    fChain->SetBranchAddress("truth_barcode"              , &truth_barcode);
    fChain->SetBranchAddress("truth_qual"                 , &truth_qual);
    fChain->SetBranchAddress("truth_pt"                 , &truth_pt);
    fChain->SetBranchAddress("truth_eta"                 , &truth_eta);
    fChain->SetBranchAddress("truth_phi"                 , &truth_phi);
    fChain->SetBranchAddress("truth_m"                 , &truth_m);
    fChain->SetBranchAddress("truth_parents"              , &truth_parents);
    fChain->SetBranchAddress("truth_children"              , &truth_children);
    fChain->SetBranchAddress("EventWeights"               , &EventWeights);
    fChain->SetBranchAddress("Q"               , &Q);

    fChain->SetBranchAddress("truth_mupair_pt1"           , &muon_pair_muon1_pt);
    fChain->SetBranchAddress("truth_mupair_eta1"          , &muon_pair_muon1_eta);
    fChain->SetBranchAddress("truth_mupair_phi1"          , &muon_pair_muon1_phi);
    fChain->SetBranchAddress("truth_mupair_ch1"           , &muon_pair_muon1_ch);
    fChain->SetBranchAddress("truth_mupair_bar1"          , &muon_pair_muon1_bar);

    fChain->SetBranchAddress("truth_mupair_pt2"           , &muon_pair_muon2_pt);
    fChain->SetBranchAddress("truth_mupair_eta2"          , &muon_pair_muon2_eta);
    fChain->SetBranchAddress("truth_mupair_phi2"          , &muon_pair_muon2_phi);
    fChain->SetBranchAddress("truth_mupair_ch2"           , &muon_pair_muon2_ch);
    fChain->SetBranchAddress("truth_mupair_bar2"          , &muon_pair_muon2_bar);

    fChain->SetBranchAddress("truth_mupair_asym"          , &truth_mupair_asym);
    fChain->SetBranchAddress("truth_mupair_acop"          , &truth_mupair_acop);
    fChain->SetBranchAddress("truth_mupair_pt"            , &truth_mupair_pt);
    fChain->SetBranchAddress("truth_mupair_y"             , &truth_mupair_y);
    // fChain->SetBranchAddress("truth_mupair_phi"           , &truth_mupair_phi);
    fChain->SetBranchAddress("truth_mupair_m"             , &truth_mupair_m);


    //SetBranch Status
    fChain->SetBranchStatus("*"                               ,0);//switch off all branches, then enable just the ones that we need
  
    fChain->SetBranchStatus("truth_id"           ,1);
    fChain->SetBranchStatus("truth_barcode"              ,1);
    fChain->SetBranchStatus("truth_qual"             ,1);
    fChain->SetBranchStatus("truth_pt"             ,1);
    fChain->SetBranchStatus("truth_eta"             ,1);
    fChain->SetBranchStatus("truth_phi"             ,1);
    fChain->SetBranchStatus("truth_m"             ,1);
    fChain->SetBranchStatus("truth_parents"             ,1);
    fChain->SetBranchStatus("truth_children"             ,1);
    fChain->SetBranchStatus("EventWeights"             ,1);
    fChain->SetBranchStatus("Q"             ,1);

    fChain->SetBranchStatus("truth_mupair_pt1"           ,1);
    fChain->SetBranchStatus("truth_mupair_eta1"              ,1);
    fChain->SetBranchStatus("truth_mupair_phi1"             ,1);
    fChain->SetBranchStatus("truth_mupair_ch1"             ,1);
    fChain->SetBranchStatus("truth_mupair_bar1"             ,1);

    fChain->SetBranchStatus("truth_mupair_pt2"             ,1);
    fChain->SetBranchStatus("truth_mupair_eta2"         ,1);
    fChain->SetBranchStatus("truth_mupair_phi2"              ,1);
    fChain->SetBranchStatus("truth_mupair_ch2"              ,1);
    fChain->SetBranchStatus("truth_mupair_bar2"              ,1);

    fChain->SetBranchStatus("truth_mupair_asym"              ,1);
    fChain->SetBranchStatus("truth_mupair_acop"           ,1);
    fChain->SetBranchStatus("truth_mupair_pt"             ,1);
    fChain->SetBranchStatus("truth_mupair_y"             ,1);
    // fChain->SetBranchStatus("truth_mupair_phi"         ,1);
    fChain->SetBranchStatus("truth_mupair_m"              ,1);   
}


void PowhegNTupleFirstPass::InitOutput(){
    InitOutputTrees();
    InitOutputHists();
}

void PowhegNTupleFirstPass::InitOutputTrees(){
    m_outfile=new TFile(Form("%s%smuon_pairs_powheg%s_%s_part%d.root", mcdir.c_str(), data_subdir.c_str(), dt_suffix.c_str(), mc_mode.c_str(), file_batch),"recreate");
    
    meta_tree = new TTree("meta_tree","meta_tree");
    meta_tree->Branch("nentries_before_cuts", &nentries, "nentries_before_cuts/L");

    if (!output_single_muon_tree){        
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            muonPairOutTree[ksign] = new TTree(Form("muon_pair_tree_sign%u",ksign+1),Form("all muon pairs, sign%u",ksign+1));
            muonPairOutTree[ksign]->Branch("MuonPairObj",&mpair_raw_ptr);
        }
    }
}

void PowhegNTupleFirstPass::InitOutputHists(){
    m_outHistFile=new TFile(Form("%s%shists_powheg_ntuple_processing_powheg%s_%s_part%d.root", mcdir.c_str(), data_subdir.c_str(), dt_suffix.c_str(), mc_mode.c_str(), file_batch),"recreate");
    DimuonAnalysisBaseClass::InitOutput(); // initialize histograms for cut acceptance
}


void PowhegNTupleFirstPass::Finalize(){
    m_outfile->Write();
    m_outHistFile->Write();
}

void PowhegNTupleFirstPass::FillMuonPair(int pair_ind, std::shared_ptr<MuonPairPowheg> const& mpair){
  // // use ind to record barcode instead
  mpair->m1.ind     = muon_pair_muon1_bar->at(pair_ind);
  mpair->m2.ind     = muon_pair_muon2_bar->at(pair_ind);
  mpair->m1.charge  = muon_pair_muon1_ch->at(pair_ind);//sign of pt stores charge
  mpair->m2.charge  = muon_pair_muon2_ch->at(pair_ind);//sign of pt stores charge

  // if charge unequal gives 1 (opposite sign); otherwise gives 0 (same sign)
  h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(double(nocut) + 0.5, mpair->weight);

  mpair->m1.pt    = fabs(muon_pair_muon1_pt->at(pair_ind))/1000.0;//pt of the first muon in the pair
  mpair->m2.pt    = fabs(muon_pair_muon2_pt->at(pair_ind))/1000.0;//pt of the second muon in the pair
  mpair->m1.eta   = muon_pair_muon1_eta->at(pair_ind);
  mpair->m2.eta   = muon_pair_muon2_eta->at(pair_ind);
  mpair->m1.phi   = muon_pair_muon1_phi->at(pair_ind);
  mpair->m2.phi   = muon_pair_muon2_phi->at(pair_ind);

  mpair->Q          = Q;
  // mpair->weight     = static_cast<double>(EventWeights->at(0) * filter_effcy / nentries);
  mpair->weight     = static_cast<double>(EventWeights->at(0) * filter_effcy);
  mpair->crossx     = EventWeights->at(0);
  mpair->minv       = truth_mupair_m->at(pair_ind)/1000.;
  mpair->pair_pt    = truth_mupair_pt->at(pair_ind)/1000.;
  mpair->pair_y     = truth_mupair_y->at(pair_ind);
  mpair->asym       = abs(truth_mupair_asym->at(pair_ind));
  mpair->acop       = abs(truth_mupair_acop->at(pair_ind));
  assert (mpair->asym >= 0 && mpair->acop >= 0);
}

bool PowhegNTupleFirstPass::PassCuts(const std::shared_ptr<MuonPair>& mpair){
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
 
  return true;
}


void PowhegNTupleFirstPass::ProcessData(){
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

            mpair = std::make_shared<MuonPairPowheg>(MuonPairPowheg());

            FillMuonPair(pair_ind, mpair);
            mpair->m1.ev_num = jentry;
            mpair->m2.ev_num = jentry;

            // ------------------------------------------------------------

            // Apply cuts

            // Trigger match for muon pair
            // if(dimuon_b_HLT_2mu4->at(pair_ind)==false) continue;

            if (abs(mpair->weight) > crossx_cut * filter_effcy) continue;
            if (!PassCuts(mpair))continue;
        
            //------------------------------------------------------------

            //Two things at this step: 
            //1) sort pt, eta, phi by pt
            //2) update the muon-pair values
            // mpair->Update();
            mpair->UpdateShort(); // only afterwards can we use mpair->same_sign

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
                        tempmuon = &(mpair->m1);
                        FillSingleMuonTree();
                  }
                  it = std::find(muon_index_list.begin(),muon_index_list.end(),mpair->m2.ind);
                  if(it == muon_index_list.end()){ //muon1 index NOT found
                        muon_index_list.push_back(mpair->m2.ind);
                        tempmuon = &(mpair->m2);
                        FillSingleMuonTree();
                  }
                }
                else{ // fill muon pair trees
                    FillMuonPairTree();
                }
            }

            if (!is_fullsim) PerformTruthPairAnalysis();

        } // finish second loop over muon pairs

        muon_pair_list_cur_event_pre_resonance_cut.clear();
    }//loop over events

}

void PowhegNTupleFirstPass::FillMuonPairTree(){

    // NECESSARY step: ALWAYS update the raw pointer with the current content of the shared pointer BEFORE FILLING OUTPUT TREES
    try{
        if (!mpair) throw std::runtime_error("FillMuonPairTree: Muon Pair doesn't exist!");
        mpair_raw_ptr = mpair.get();
    }catch(const std::runtime_error& e){
        std::cout << "Runtime error caught in function FillMuonPairTree: " << e.what() << std::endl;
        std::cout << "Return without filling the muon pair in the output trees!" << std::endl;
        return;
    }

    int nsign = (mpair->same_sign)? 0:1;
    muonPairOutTree[nsign]->Fill();
}


void PowhegNTupleFirstPass::HistAdjust(){
  DimuonAnalysisBaseClass::HistAdjust();
}




#endif
