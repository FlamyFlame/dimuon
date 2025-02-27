#ifndef __PbPbDataNTupleFirstPass_C__
#define __PbPbDataNTupleFirstPass_C__

#include "PbPbDataNTupleFirstPass.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"

//initialize the TChain
void PbPbDataNTupleFirstPass::InitInput(){

   fChain = new TChain("HeavyIonD3PD","HeavyIonD3PD");
   fChain->SetMakeClass(1);
   fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2/All_Data2018_12March2022.root");
   fChain->Add("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2/All_Data2015_12March2022.root"); 

   fChain->SetBranchAddress("RunNumber"                   , &RunNumber);
   fChain->SetBranchAddress("ZdcEtA"                      , &ZdcEtA              );
   fChain->SetBranchAddress("ZdcEtC"                      , &ZdcEtC              );
   fChain->SetBranchAddress("FCal_Et"                     , &FCal_Et);
   fChain->SetBranchAddress("centrality"                  , &centrality);
   fChain->SetBranchAddress("muon_deltaP_overP"           , &muon_deltaP_overP);

   fChain->SetBranchAddress("muon_pair_muon1_index"       , &muon_pair_muon1_index);
   fChain->SetBranchAddress("muon_pair_muon1_pt"          , &muon_pair_muon1_pt);
   fChain->SetBranchAddress("muon_pair_muon1_eta"         , &muon_pair_muon1_eta);
   fChain->SetBranchAddress("muon_pair_muon1_phi"         , &muon_pair_muon1_phi);
   fChain->SetBranchAddress("muon_pair_muon1_quality"     , &muon_pair_muon1_quality);
   fChain->SetBranchAddress("muon_pair_muon1_d0"          , &muon_pair_muon1_d0);
   fChain->SetBranchAddress("muon_pair_muon1_z0"          , &muon_pair_muon1_z0);

   fChain->SetBranchAddress("muon_pair_muon2_index"       , &muon_pair_muon2_index);
   fChain->SetBranchAddress("muon_pair_muon2_pt"          , &muon_pair_muon2_pt);
   fChain->SetBranchAddress("muon_pair_muon2_eta"         , &muon_pair_muon2_eta);
   fChain->SetBranchAddress("muon_pair_muon2_phi"         , &muon_pair_muon2_phi);
   fChain->SetBranchAddress("muon_pair_muon2_quality"     , &muon_pair_muon2_quality);
   fChain->SetBranchAddress("muon_pair_muon2_d0"          , &muon_pair_muon2_d0);
   fChain->SetBranchAddress("muon_pair_muon2_z0"          , &muon_pair_muon2_z0);
   fChain->SetBranchAddress("b_HLT_mu4_mu4noL1"           , &b_HLT_mu4_mu4noL1);
   fChain->SetBranchAddress("dimuon_b_HLT_mu4_mu4noL1"    , &dimuon_b_HLT_mu4_mu4noL1);


   //SetBranch Status
   fChain->SetBranchStatus("*"                               ,0);//switch off all branches, then enable just the ones that we need
   fChain->SetBranchStatus("RunNumber"                       ,1);
   fChain->SetBranchStatus("ZdcEtA"                          ,1);
   fChain->SetBranchStatus("ZdcEtC"                          ,1);
   fChain->SetBranchStatus("FCal_Et"                         ,1);
   fChain->SetBranchStatus("centrality"                      ,1);
   fChain->SetBranchStatus("muon_deltaP_overP"               ,1);

   fChain->SetBranchStatus("muon_pair_muon1_index"           ,1);
   fChain->SetBranchStatus("muon_pair_muon1_pt"              ,1);
   fChain->SetBranchStatus("muon_pair_muon1_eta"             ,1);
   fChain->SetBranchStatus("muon_pair_muon1_phi"             ,1);
   fChain->SetBranchStatus("muon_pair_muon1_quality"         ,1);
   fChain->SetBranchStatus("muon_pair_muon1_d0"              ,1);
   fChain->SetBranchStatus("muon_pair_muon1_z0"              ,1);

   fChain->SetBranchStatus("muon_pair_muon2_index"           ,1);
   fChain->SetBranchStatus("muon_pair_muon2_pt"              ,1);
   fChain->SetBranchStatus("muon_pair_muon2_eta"             ,1);
   fChain->SetBranchStatus("muon_pair_muon2_phi"             ,1);
   fChain->SetBranchStatus("muon_pair_muon2_quality"         ,1);
   fChain->SetBranchStatus("muon_pair_muon2_d0"              ,1);
   fChain->SetBranchStatus("muon_pair_muon2_z0"              ,1);

   fChain->SetBranchStatus("b_HLT_mu4_mu4noL1"               ,1);
   fChain->SetBranchStatus("dimuon_b_HLT_mu4_mu4noL1"        ,1);
}

void PbPbDataNTupleFirstPass::InitOutput(){

  if (mode == 3){
    m_outfile=new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2/single_muon_trees_small_ctr_intvls_pbpb_run2.root","recreate");
    muonOutTree = new TTree("muon_tree","all single muons");
    muonOutTree->Branch("MuonObj",&tempmuon);

    for (int jctr = 0; jctr < ParamsSet::nCtrIntvls; jctr++){
      muonOutTreeBinned[jctr] = new TTree(Form("muon_tree_ctr%d",jctr+1),Form("all muons, centrality bin %d",jctr+1));
      muonOutTreeBinned[jctr]->Branch("MuonObj",&tempmuon);
    }
  }
  else{ // mode == 4
    m_outfile=new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2/muon_pairs_small_ctr_intvls_pbpb_run2.root","recreate");

    for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
      muonPairOutTree[ksign] = new TTree(Form("muon_pair_tree_sign%u",ksign+1),Form("all muon pairs, sign%u",ksign+1));
      muonPairOutTree[ksign]->Branch("MuonPairObj",&mpair_raw_ptr);
      for (unsigned int jctr = 0; jctr < ParamsSet::nCtrIntvls; jctr++){
        muonPairOutTreeBinned[jctr][ksign] = new TTree(Form("muon_pair_tree_ctr%u_sign%u",jctr+1,ksign+1),Form("all muon pairs, ctr%u, sign%u",jctr+1,ksign+1));
        muonPairOutTreeBinned[jctr][ksign]->Branch("MuonPairObj",&mpair_raw_ptr);
      }
    }
  }
  
  // ---------------------------------------------------------------------------------------------------------------------------

  m_outHistFile=new TFile("/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_run2/hists_cut_acceptance_pbpb_run2.root","recreate");
  MuonNTupleFirstPassBaseClass::InitOutput();
}

void PbPbDataNTupleFirstPass::FillMuonPair(int pair_ind, std::shared_ptr<MuonPairData>& mpair){
  mpair->m1.ind     = muon_pair_muon1_index->at(pair_ind);
  mpair->m2.ind     = muon_pair_muon2_index->at(pair_ind);

  mpair->m1.pt    = fabs(muon_pair_muon1_pt->at(pair_ind))/1000.0;//pt of the first muon in the pair
  mpair->m2.pt    = fabs(muon_pair_muon2_pt->at(pair_ind))/1000.0;//pt of the second muon in the pair
  mpair->m1.eta   = muon_pair_muon1_eta->at(pair_ind);
  mpair->m2.eta   = muon_pair_muon2_eta->at(pair_ind);
  mpair->m1.phi   = muon_pair_muon1_phi->at(pair_ind);
  mpair->m2.phi   = muon_pair_muon2_phi->at(pair_ind);
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
  mpair->m1.ev_centrality = centrality;
  mpair->m2.ev_centrality = centrality;
  mpair->m1.ev_FCal_Et = FCal_Et;
  mpair->m2.ev_FCal_Et = FCal_Et;
}

bool PbPbDataNTupleFirstPass::PassCuts(const std::shared_ptr<MuonPair>& mpair){
  //Apply ALL CUTS but for resonances

  //require some quality cuts on the muons
  if((mpair->m1.quality&mpair->m2.quality&1  )==0) return false;//combined muon
  if((mpair->m1.quality&mpair->m2.quality&8  )==0) return false;//Medium muon
  if((mpair->m1.quality&mpair->m2.quality&32 )==0) return false;//IDCuts
  if((mpair->m1.quality&mpair->m2.quality&256)==0) return false;//MuonCuts
  h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(pass_muon_quality + 0.5, mpair->weight); // if same sign: fill the h_cutAcceptance[0] histogram; if opposite sign, fill the h_cutAcceptance[1] histogram

  //cut on z0pair - order doesn't mpair->m1tter
  // double z0sinTheta1 = mpair->m1.z0*sin(2.0*atan(exp(-mpair->m1.eta)));
  // double z0sinTheta2 = mpair->m2.z0*sin(2.0*atan(exp(-mpair->m2.eta)));
  // float  z0_combined =sqrt(z0sinTheta1*z0sinTheta1 + z0sinTheta2*z0sinTheta2);
  // if(z0_combined>1.0) return false;

  if (fabs(mpair->m1.eta) > 2.4 || fabs(mpair->m2.eta) > 2.4) return false;
  h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(pass_muon_eta + 0.5, mpair->weight);

  if (mpair->m1.pt < 4 || mpair->m2.pt < 4) return false;
  h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(pass_muon_pt + 0.5, mpair->weight);
  
  if( fabs(mpair->m1.dP_overP) > pms.deltaP_overP_thrsh || fabs(mpair->m2.dP_overP) > pms.deltaP_overP_thrsh ) return false;
  h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(pass_muon_dP_overP + 0.5, mpair->weight);
 
  return true;
}

void PbPbDataNTupleFirstPass::FillSingleMuonTree(){
  muonOutTree->Fill();
  for (int jctr = 0; jctr < pms.nCtrIntvls; jctr++){
    if (tempmuon->ev_centrality >= pms.CtrStep * jctr && tempmuon->ev_centrality < pms.CtrStep *(jctr+1)){
      muonOutTreeBinned[jctr]->Fill();
    }
  }
}

void PbPbDataNTupleFirstPass::FillMuonPairTree(){

  // NECESSARY step: ALWAYS update the raw pointer with the current content of the shared pointer BEFORE FILLING OUTPUT TREES
  try{
    if (!mpair) throw std::runtime_error();
    mpair_raw_ptr = mpair.get();
  }catch(const std::runtime_error& e){
    std::cout << "Runtime error caught in function FillMuonPairTree: " << e.what() << std::endl;
    std::cout << "Return without filling the muon pair in the output trees!" << std::endl;
    return;
  }

  int nsign = (mpair->same_sign)? 0:1;
    muonPairOutTree[nsign]->Fill();
    for (unsigned int jctr = 0; jctr < pms.nCtrIntvls; jctr++){
      if (mpair->avg_centrality >= pms.CtrStep * jctr && mpair->avg_centrality < pms.CtrStep *(jctr+1)){
        muonPairOutTreeBinned[jctr][nsign]->Fill();
      }
    }
}

void PbPbDataNTupleFirstPass::ProcessData(){

  Long64_t nentries = fChain->GetEntries();//number of events
  // Long64_t nentries = 100000;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {//loop over the events

    if(jentry%100000==0) cout<<"Processing "<<jentry<<" event out of "<<nentries<<" events"<<std::endl;

    int num_bytes = fChain->GetEntry(jentry);//read in an event
    if(num_bytes==0){
      std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
      throw std::exception();
    }
 
    muon_pair_list_cur_event_pre_resonance_cut.clear();
    resonance_tagged_muon_index_list.clear(); // MUST CLEAR for each event!!

    //trigger requirement for event
    if(!b_HLT_mu4_mu4noL1) continue;


    std::vector<int> muon_index_list = {};
    std::vector<int>::iterator it;

    int NPairs=muon_pair_muon1_pt->size();//number of muon pairs in the event

    for(int pair_ind=0;pair_ind<NPairs;pair_ind++){//first loop over all muon-pairs in the event

      mpair = std::make_shared<MuonPairData>(MuonPairData());

      FillMuonPair(pair_ind, mpair);
      mpair->m1.ev_num = jentry;
      mpair->m2.ev_num = jentry;

      //------------------------------------------------------------

      //Apply cuts
      
      h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(nocut + 0.5, mpair->weight);

      //Trigger match for muon pair
      if(dimuon_b_HLT_mu4_mu4noL1->at(pair_ind)==false) continue;

      if (!PassCuts(mpair))continue;
    
      //------------------------------------------------------------

      //Two things at this step: 
      //1) sort pt, eta, phi by pt
      //2) update the muon-pair values
      mpair->Update();

      // resonance tag
      ResonanceTagging(mpair);

      // photo-production cut
      if (IsPhotoProduction(mpair)) continue;
      h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(pass_photoprod + 0.5, mpair->weight);

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
      else{// mode = 4: fill muon pair trees
        FillMuonPairTree();
      }
    }
  }//loop over events
}

void PbPbDataNTupleFirstPass::Run(){

  if(mode != 3 && mode != 4){
    std::cout<<"Error:: Mode has to equal 3 or 4; code is used for outputting muon / muon-pair trees only."<<std::endl;
    throw std::exception();
  }

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  InitInput();
  InitOutput();
  ProcessData();
  HistAdjust();

  m_outfile->Write();
  m_outHistFile->Write();

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "CPU time used (in seconds) for mode " << mode << " is " << cpu_time_used << std::endl;

}

#endif
