#include "PPDataNTupleFirstPass.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"

PPDataNTupleFirstPass::PPDataNTupleFirstPass(){
  numCuts = numCuts_data;
  cutLabels = cutLabels_data;

  isTight = false;
  isRun3 = true;
  run3_file_batch = 0;
  run3_use_mu4_mu4_noL1 = true;
  
  std::cout << "PP Data Ntuple processing script:" << std::endl;
  std::cout << "The following public variable(s) **MUST** be set:" << std::endl;
  std::cout << "--> run3_file_batch: int that decides which run3-file batch to process, only has effect when isRun3 is true" << std::endl;
  std::cout << "                     value = 1,2,3,4 (DEFAULT 0 --> MUST BE SET)" << std::endl;
  std::cout << std::endl;

  std::cout << "The following public variable(s) should be checked:" << std::endl;
  std::cout << "--> resonance_cut_mode: integer that determines which set of resonant cuts to apply" << std::endl;
  std::cout << "        * resonance_cut_mode = 0: NO resonance cut" << std::endl;
  std::cout << "        * resonance_cut_mode = 1: old resonance cut" << std::endl;
  std::cout << "        * resonance_cut_mode = 1: new resonance cut (default)" << std::endl;
  std::cout << "        If resonance_cut_mode value is outside {0,1,2}: assume default option" << std::endl;
  std::cout << "--> isTight: boolean, default false - if true: require tight WP; false; require medium WP" << std::endl;
  std::cout << "--> isRun3: boolean, default true - if true: run run3 data; false: run run2 data" << std::endl;
  std::cout << "--> run3_use_mu4_mu4_noL1: boolean, default true, only has effect when isRun3 is true" << std::endl;
  std::cout << "                           if true: requires mu4_mu4_noL1 trigger; false: requires 2mu4 trigger" << std::endl;
  std::cout << std::endl;

  std::cout << "if isRun3, output files will be written to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024" << std::endl;
  std::cout << "else,      output files will be written to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_run2" << std::endl;
  std::cout << "" << std::endl;
}

//initialize the TChain
void PPDataNTupleFirstPass::InitInput(){

    fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
    fChain->SetMakeClass(1);
    if (!isRun3){ //run2
        fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/New_All_DataPP2017_5TeV_Dec2021.root");
    }else{ //run3
        switch (run3_file_batch){
        case 1:
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part1.root");
            break;
        case 2:
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part3/data_pp24_part3_1.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part3/data_pp24_part3_2.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part3/data_pp24_part3_3.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part3/data_pp24_part3_4.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part3/data_pp24_part3_5.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part3/data_pp24_part3_6.root");
            break;
        case 3:
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_1.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_2.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_4.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_5.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_6.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_8.root");
            break;
        case 4:
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_3.root");
            fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/data_pp24_part2/data_pp24_part2_7.root");
            break;
        default:
            throw std::runtime_error("Run3 file batch invalid/unspecified: have to be between 1 and 4. No result gets run.");
        }
    }
  
    fChain->SetBranchAddress("RunNumber"                   , &RunNumber);
    fChain->SetBranchAddress("muon_deltaP_overP"           , &muon_deltaP_overP);
  
    fChain->SetBranchAddress("muon_pair_muon1_index"       , &muon_pair_muon1_index);
    fChain->SetBranchAddress("muon_pair_muon1_pt"          , &muon_pair_muon1_pt);
    fChain->SetBranchAddress("muon_pair_muon1_eta"         , &muon_pair_muon1_eta);
    fChain->SetBranchAddress("muon_pair_muon1_phi"         , &muon_pair_muon1_phi);
    fChain->SetBranchAddress("muon_pair_muon1_quality"     , &muon_pair_muon1_quality);
    fChain->SetBranchAddress("muon_pair_muon1_d0"          , &muon_pair_muon1_d0);
    fChain->SetBranchAddress("muon_pair_muon1_z0"          , &muon_pair_muon1_z0);
    fChain->SetBranchAddress("muon_pair_muon1_trk_pt"      , &muon_pair_muon1_trk_pt);
    fChain->SetBranchAddress("muon_pair_muon1_trk_eta"     , &muon_pair_muon1_trk_eta);
    fChain->SetBranchAddress("muon_pair_muon1_trk_phi"     , &muon_pair_muon1_trk_phi);
  
    fChain->SetBranchAddress("muon_pair_muon2_index"       , &muon_pair_muon2_index);
    fChain->SetBranchAddress("muon_pair_muon2_pt"          , &muon_pair_muon2_pt);
    fChain->SetBranchAddress("muon_pair_muon2_eta"         , &muon_pair_muon2_eta);
    fChain->SetBranchAddress("muon_pair_muon2_phi"         , &muon_pair_muon2_phi);
    fChain->SetBranchAddress("muon_pair_muon2_quality"     , &muon_pair_muon2_quality);
    fChain->SetBranchAddress("muon_pair_muon2_d0"          , &muon_pair_muon2_d0);
    fChain->SetBranchAddress("muon_pair_muon2_z0"          , &muon_pair_muon2_z0);
    fChain->SetBranchAddress("muon_pair_muon2_trk_pt"      , &muon_pair_muon2_trk_pt);
    fChain->SetBranchAddress("muon_pair_muon2_trk_eta"     , &muon_pair_muon2_trk_eta);
    fChain->SetBranchAddress("muon_pair_muon2_trk_phi"     , &muon_pair_muon2_trk_phi);

    if (!isRun3){ // run2
        // no mu4_mu4_noL1 trigger in run2
        // dimuon_b_HLT_2mu4 is a vector of boolean containing trigger match for each muon pair
        fChain->SetBranchAddress("b_HLT_2mu4"                  , &b_HLT_2mu4);
        fChain->SetBranchAddress("dimuon_b_HLT_2mu4"           , &dimuon_b_HLT_2mu4);        
    } else{ // run3
        // we need both the 2mu4 trigger (requiring 2 muons at L1 and at HLT)
        // and the mu4_mu4_noL1 trigger (requiring 1 muon at L1 and 2 muons at HLT)
        // since (1) L1MU3V is prescaled at L1 (2) mu4_mu4_noL1 as a supporting trigger is further prescaled at HLT
        fChain->SetBranchAddress("b_HLT_2mu4_L12MU3V"                       , &b_HLT_2mu4);
        fChain->SetBranchAddress("b_HLT_mu4_mu4noL1_L1MU3V"                 , &b_HLT_mu4_mu4noL1);
        fChain->SetBranchAddress("dimuon_b_HLT_2mu4_L12MU3V"                , &dimuon_b_HLT_2mu4);
        fChain->SetBranchAddress("dimuon_b_HLT_mu4_mu4noL1_L1MU3V"          , &dimuon_b_HLT_mu4_mu4noL1);
    }

    //SetBranch Status
    fChain->SetBranchStatus("*"                               ,0);//switch off all branches, then enable just the ones that we need
    fChain->SetBranchStatus("RunNumber"                       ,1);
    fChain->SetBranchStatus("muon_deltaP_overP"               ,1);
  
    fChain->SetBranchStatus("muon_pair_muon1_index"           ,1);
    fChain->SetBranchStatus("muon_pair_muon1_pt"              ,1);
    fChain->SetBranchStatus("muon_pair_muon1_eta"             ,1);
    fChain->SetBranchStatus("muon_pair_muon1_phi"             ,1);
    fChain->SetBranchStatus("muon_pair_muon1_trk_pt"              ,1);
    fChain->SetBranchStatus("muon_pair_muon1_trk_eta"             ,1);
    fChain->SetBranchStatus("muon_pair_muon1_trk_phi"             ,1);
    fChain->SetBranchStatus("muon_pair_muon1_quality"         ,1);
    fChain->SetBranchStatus("muon_pair_muon1_d0"              ,1);
    fChain->SetBranchStatus("muon_pair_muon1_z0"              ,1);
  
    fChain->SetBranchStatus("muon_pair_muon2_index"           ,1);
    fChain->SetBranchStatus("muon_pair_muon2_pt"              ,1);
    fChain->SetBranchStatus("muon_pair_muon2_eta"             ,1);
    fChain->SetBranchStatus("muon_pair_muon2_phi"             ,1);
    fChain->SetBranchStatus("muon_pair_muon2_trk_pt"              ,1);
    fChain->SetBranchStatus("muon_pair_muon2_trk_eta"             ,1);
    fChain->SetBranchStatus("muon_pair_muon2_trk_phi"             ,1);
    fChain->SetBranchStatus("muon_pair_muon2_quality"         ,1);
    fChain->SetBranchStatus("muon_pair_muon2_d0"              ,1);
    fChain->SetBranchStatus("muon_pair_muon2_z0"              ,1);
  
    if (!isRun3){ // run2
        fChain->SetBranchStatus("b_HLT_2mu4"               ,1);
        fChain->SetBranchStatus("dimuon_b_HLT_2mu4"        ,1);
    } else{ // run3
        fChain->SetBranchStatus("b_HLT_2mu4_L12MU3V"                        ,1);
        fChain->SetBranchStatus("b_HLT_mu4_mu4noL1_L1MU3V"                  ,1);
        fChain->SetBranchStatus("dimuon_b_HLT_2mu4_L12MU3V"                 ,1);
        fChain->SetBranchStatus("dimuon_b_HLT_mu4_mu4noL1_L1MU3V"           ,1);
    }
}

void PPDataNTupleFirstPass::InitOutput(){

    std::string tight_suffix = (isTight)? "_tight" : "";
    
    std::string output_dir = "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/";
    std::string run_dir = isRun3? "pp_2024/" : "pp_run2/";
    output_dir += run_dir;

    std::string run_suffix = isRun3? "_pp_2024" : "_pp_run2";

    std::map<int, std::string> mode_to_file_name_map;
    mode_to_file_name_map[3] = "single_muon_trees";
    mode_to_file_name_map[4] = "muon_pairs";

    std::string run3_trig_suffix = "";
    if (isRun3){
      run3_trig_suffix = run3_use_mu4_mu4_noL1? "_mu4_mu4noL1" : "_2mu4";
    }

    std::map<int, std::string> run3_batch_suffix_map;
    run3_batch_suffix_map[1] = "_part1";
    run3_batch_suffix_map[2] = "_part3";
    run3_batch_suffix_map[3] = "_part2-1";
    run3_batch_suffix_map[4] = "_part2-2";

    std::string run3_batch_suffix = "";
    if (isRun3) run3_batch_suffix = run3_batch_suffix_map[run3_file_batch];

    std::string resonance_cut_suffix = "";
    switch (resonance_cut_mode){
    case 0:
      resonance_cut_suffix = "_no_res_cut";
      break;
    case 1:
      resonance_cut_suffix = "_old_res_cut";
      break;
    case 2:
      resonance_cut_suffix = "_new_res_cut";
      break;
    default:
      std::cout << "Public variable resonance_cut_mode is set to a value outside {0,1,2}: INVALID. Apply new resonance cuts by default." << std::endl;
      resonance_cut_suffix = "_new_res_cut";
    }

    output_file_path = output_dir + mode_to_file_name_map[mode] + run_suffix + run3_batch_suffix + run3_trig_suffix + tight_suffix + resonance_cut_suffix + ".root";
    output_hist_file_path = output_dir + "hists_cut_acceptance" + run_suffix + run3_batch_suffix + run3_trig_suffix + tight_suffix + resonance_cut_suffix + ".root";

    // ------------------------------------------------------------------------------

    m_outfile=new TFile(output_file_path.c_str(),"recreate");
    
    if (mode == 3){
        muonOutTree = new TTree("muon_tree","all single muons");
        muonOutTree->Branch("MuonObj",&tempmuon);
    }
    else{ // mode == 4: output muon pair trees
        for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
            muonPairOutTree[ksign] = new TTree(Form("muon_pair_tree_sign%u",ksign+1),Form("all muon pairs, sign%u",ksign+1));
            muonPairOutTree[ksign]->Branch("MuonPairObj",&mpair_raw_ptr);
        }
    }

    // ------------------------------------------------------------------------------

    m_outHistFile=new TFile(output_hist_file_path.c_str(),"recreate");
    MuonNTupleFirstPassBaseClass::InitOutput();
}

void PPDataNTupleFirstPass::FillMuonPair(int pair_ind, std::shared_ptr<MuonPairData>& mpair){
  mpair->m1.ind     = muon_pair_muon1_index->at(pair_ind);
  mpair->m2.ind     = muon_pair_muon2_index->at(pair_ind);

  mpair->m1.pt    = fabs(muon_pair_muon1_pt->at(pair_ind))/1000.0;//pt of the first muon in the pair
  mpair->m2.pt    = fabs(muon_pair_muon2_pt->at(pair_ind))/1000.0;//pt of the second muon in the pair
  mpair->m1.eta   = muon_pair_muon1_eta->at(pair_ind);
  mpair->m2.eta   = muon_pair_muon2_eta->at(pair_ind);
  mpair->m1.phi   = muon_pair_muon1_phi->at(pair_ind);
  mpair->m2.phi   = muon_pair_muon2_phi->at(pair_ind);

  mpair->m1_trk_pt    = fabs(muon_pair_muon1_trk_pt->at(pair_ind))/1000.0;//pt of the first muon in the pair
  mpair->m2_trk_pt    = fabs(muon_pair_muon2_trk_pt->at(pair_ind))/1000.0;//pt of the second muon in the pair
  mpair->m1_trk_eta   = muon_pair_muon1_trk_eta->at(pair_ind);
  mpair->m2_trk_eta   = muon_pair_muon2_trk_eta->at(pair_ind);
  mpair->m1_trk_phi   = muon_pair_muon1_trk_phi->at(pair_ind);
  mpair->m2_trk_phi   = muon_pair_muon2_trk_phi->at(pair_ind);

  mpair->m1.d0    = muon_pair_muon1_d0 ->at(pair_ind);
  mpair->m2.d0    = muon_pair_muon2_d0 ->at(pair_ind);
  mpair->m1.z0    = muon_pair_muon1_z0 ->at(pair_ind);
  mpair->m2.z0    = muon_pair_muon2_z0 ->at(pair_ind);
  mpair->m1.charge  =(muon_pair_muon1_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge
  mpair->m2.charge  =(muon_pair_muon2_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge
  mpair->m1.quality = muon_pair_muon1_quality->at(pair_ind);
  mpair->m2.quality = muon_pair_muon2_quality->at(pair_ind);
  mpair->m1.dP_overP = muon_deltaP_overP->at(mpair->m1.ind);
  mpair->m2.dP_overP = muon_deltaP_overP->at(mpair->m2.ind);
}

bool PPDataNTupleFirstPass::PassCuts(const std::shared_ptr<MuonPair>& mpair){
  //Apply ALL CUTS but for resonances

  //require some quality cuts on the muons
  if((mpair->m1.quality&mpair->m2.quality&1  )==0) return false;//compair->m2ined muon
  if (isTight){
    if ((mpair->m1.quality&mpair->m2.quality&16  )==0) return false;//Tight muon
  }else{
    if ((mpair->m1.quality&mpair->m2.quality&8  )==0) return false;//Medium muon
  }
  if((mpair->m1.quality&mpair->m2.quality&32 )==0) return false;//IDCuts
  //if((mpair->m1.quality&mpair->m2.quality&256)==0) return false;//MuonCuts
  
  // if pass quality cut, fill in the pair with weight in the cut-acceptance histogram
  h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(pass_muon_quality + 0.5, mpair->weight);

  if (fabs(mpair->m1.eta) > 2.4 || fabs(mpair->m2.eta) > 2.4) return false;
  h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(pass_muon_eta + 0.5, mpair->weight);

  if (mpair->m1.pt < 4 || mpair->m2.pt < 4) return false;
  h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(pass_muon_pt + 0.5, mpair->weight);

  if (fabs(mpair->m1.dP_overP) > pms.deltaP_overP_thrsh || fabs(mpair->m2.dP_overP) > pms.deltaP_overP_thrsh) return false;
  h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(pass_muon_dP_overP + 0.5, mpair->weight);
 
  return true;
}


bool PPDataNTupleFirstPass::TrigMatch(int pair_ind){
  //Trigger match for muon pair
  bool trig_match_req = false;
  try{
    if (isRun3 && run3_use_mu4_mu4_noL1){
      if (!dimuon_b_HLT_mu4_mu4noL1){
        throw std::runtime_error("The pointer dimuon_b_HLT_mu4_mu4noL1 (to vector of bool) is a null pointer!");
      }
      trig_match_req = dimuon_b_HLT_mu4_mu4noL1->at(pair_ind);
    }else{
      if (!dimuon_b_HLT_2mu4){
        throw std::runtime_error("The pointer dimuon_b_HLT_2mu4 (to vector of bool) is a null pointer!");
      }
      trig_match_req = dimuon_b_HLT_2mu4->at(pair_ind);
    }
  }catch (const std::out_of_range& e){
    std::cerr << "out-of-range error occured: " << e.what() << endl;
    std::cout << "By default, event will not be saved." << endl;
    trig_match_req = false;
  }catch (const std::runtime_error& e){
    std::cerr << "Run time error occured: " << e.what() << endl;
    std::cout << "By default, event will not be saved." << endl;
    trig_match_req = false;
  }catch(...){
    std::cerr << "Some other unknown error occured at the trigger matching stage" << endl;
    std::cout << "By default, event will not be saved." << endl;
    trig_match_req = false;
  }

  return trig_match_req;
}

void PPDataNTupleFirstPass::FillSingleMuonTree(){
  muonOutTree->Fill();
}

void PPDataNTupleFirstPass::FillMuonPairTree(){

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


void PPDataNTupleFirstPass::ProcessData(){

  // Long64_t nentries = 100000;
  Long64_t nentries = fChain->GetEntries();//number of events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {//loop over the events
  // for (Long64_t jentry=0; jentry<1000;jentry++) {//loop over the events

    if(jentry%100000==0) cout<<"Processing "<<jentry<<" event out of "<<nentries<<" events"<<std::endl;

    int num_bytes = fChain->GetEntry(jentry);//read in an event
    if(num_bytes==0){
      std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
      throw std::exception();
    }

    muon_pair_list_cur_event_pre_resonance_cut.clear();
    resonance_tagged_muon_index_list.clear(); // MUST CLEAR for each event!!
    resonance_tagged_muon_index_list_old.clear(); // MUST CLEAR for each event!!

    //trigger requirement for event
    bool trigger_req = (isRun3 && run3_use_mu4_mu4_noL1)? b_HLT_mu4_mu4noL1 : b_HLT_2mu4;
    if(!trigger_req) continue;

    std::vector<int> muon_index_list = {};
    std::vector<int>::iterator it;

    int NPairs=muon_pair_muon1_pt->size();//number of muon pairs in the event

    for(int pair_ind=0;pair_ind<NPairs;pair_ind++){//first loop over all muon-pairs in the event

      mpair = std::make_shared<MuonPairData>(MuonPairData());

      FillMuonPair(pair_ind, mpair);
      mpair->m1.ev_num = jentry;
      mpair->m2.ev_num = jentry;

      //------------------------------------------------------------
      h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(nocut + 0.5, mpair->weight);

      if(!TrigMatch(pair_ind)) continue;
      h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(pass_trigger_match + 0.5, mpair->weight);

      if (!PassCuts(mpair))continue;
    
      //------------------------------------------------------------

      //Two things at this step: 
      //1) sort pt, eta, phi by pt
      //2) update the muon-pair values
      // std::cout << mpair->m1.ev_num << std::endl;
      mpair->Update();

      // resonance tag
      ResonanceTagging(mpair);
      ResonanceTaggingOld(mpair);

      // photo-production cut - do NOT apply for pp
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

      // apply resonance cuts if resonance_cut_mode != 0
      if (resonance_cut_mode == 1){ // apply old cut
        itres_m1 = std::find(resonance_tagged_muon_index_list_old.begin(),resonance_tagged_muon_index_list_old.end(),mpair->m1.ind);
        if(itres_m1 != resonance_tagged_muon_index_list_old.end())  continue;

        itres_m2 = std::find(resonance_tagged_muon_index_list_old.begin(),resonance_tagged_muon_index_list_old.end(),mpair->m2.ind);
        if(itres_m2 != resonance_tagged_muon_index_list_old.end())  continue;
      } else if (resonance_cut_mode == 2){ // apply new cut
        itres_m1 = std::find(resonance_tagged_muon_index_list.begin(),resonance_tagged_muon_index_list.end(),mpair->m1.ind);
        if(itres_m1 != resonance_tagged_muon_index_list.end())  continue;

        itres_m2 = std::find(resonance_tagged_muon_index_list.begin(),resonance_tagged_muon_index_list.end(),mpair->m2.ind);
        if(itres_m2 != resonance_tagged_muon_index_list.end())  continue;
      }

      h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(pass_resonance + 0.5, mpair->weight);

      //------------------------------------------------------------

      if(mode == 3){ // single muon information; no dR selection
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
      else FillMuonPairTree(); // mode = 4
    }
  }//loop over events
}


void PPDataNTupleFirstPass::Run(){

  if(mode != 3 && mode != 4){
      std::cout<<"Error:: Mode has to equal 3 or 4; code is used for outputting muon / muon-pair trees only."<<std::endl;
      throw std::exception();
  }

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  std::cout << "Mode = " << mode << ". Output file is " << m_outfile << std::endl;
  InitInput();
  InitOutput();
  ProcessData();
  HistAdjust();

  m_outfile->Write();
  m_outHistFile->Write();
  std::cout << "Output muon-pair trees have been written to: " << output_file_path << std::endl;
  std::cout << "Output histograms have been written to: " << output_hist_file_path << std::endl;

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "CPU time used (in seconds) for mode " << mode << " is " << cpu_time_used << std::endl;

}
