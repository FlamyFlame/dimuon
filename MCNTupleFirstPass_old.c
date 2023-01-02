#ifndef __MCNTupleFirstPass_C__
#define __MCNTupleFirstPass_C__

#include "MCNTupleFirstPass.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"
#include<bits/stdc++.h>

// Benefits of using enum: like labeling equations 
// and referring to them using the labels rather than numbers
// we do not need to remember which number corresponds to which cut
// also, if we add new cuts, we don't need to change the part of the code
// corresponding to the existing cuts
enum cuts_MC{
  nocut, 
  pass_muon_eta, 
  pass_muon_pt, 
  pass_resonance, 
  pass_photoprod
};

bool MCNTupleFirstPass::PassCuts(){
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
  h_cutAcceptance[!mpair->same_sign]->Fill(pass_muon_eta + 0.5); // if same sign: fill the h_cutAcceptance[0] histogram; if opposite sign, fill the h_cutAcceptance[1] histogram

  if (mpair->m1.pt < 4 || mpair->m2.pt < 4) return false;
  
  // passed muon pt cut
  // if (mpair->same_sign) h_cutAcceptance[0]->Fill(pass_muon_pt + 0.5);
  // else h_cutAcceptance[1]->Fill(pass_muon_pt + 0.5);
  h_cutAcceptance[!mpair->same_sign]->Fill(pass_muon_pt + 0.5);

  // if( fabs(mpair->m1.dP_overP) > pms.deltaP_overP_thrsh || fabs(mpair->m2.dP_overP) > pms.deltaP_overP_thrsh ) pass = false;
 
  return true;
}

bool MCNTupleFirstPass::IsResonance(){
  // assuming we have already filled up the muon pair mpair
  if (mpair->same_sign) h_cutAcceptance[0]->Fill(pass_resonance + 0.5);
  else{ // opposite sign
    if (mpair->minv > pms.minv_upper) return true; // upper cut at 60 GeV

    // bool isresonance = false;
    
    for (array<float,2> ires : pms.minv_cuts){
      if (mpair->minv > ires[0] && mpair->minv < ires[1]){
        return true;
      }
    }
    
    h_cutAcceptance[1]->Fill(pass_resonance + 0.5);
  }
  
  return false;
}

bool MCNTupleFirstPass::IsPhotoProduction(){
  // A := |Delta PT| / (sum pT) < 0.05 && alpha := (pi-Dphi)/pi < 0.01
  // float A = (mpair->m1.pt - mpair->m2.pt) / (mpair->m1.pt + mpair->m2.pt);
  // assert (A >= 0);
  // float alpha = (pms.PI - fabs(mpair->dphi)) / pms.PI;
  // return (!(mpair->same_sign) && A < 0.05 && alpha < 0.01);
  // return (!(mpair->same_sign) && mpair->asym < 0.05 && mpair->acop < 0.01);
  
  if (mpair->same_sign) h_cutAcceptance[0]->Fill(pass_photoprod + 0.5);
  else if (mpair->asym < 0.05 && mpair->acop < 0.01) return true; // is photoproduction 
  else h_cutAcceptance[1]->Fill(pass_photoprod + 0.5);

  return false;
}

void MCNTupleFirstPass::FillSingleMuonTree(){
  muonOutTree->Fill();
}

int MCNTupleFirstPass::ParentGrouping(int parent_id){
  
  if ((parent_id >= 11 && parent_id <= 16) || parent_id == 22) // leptons, neutrinoes, photons
    return 1;
  if ((parent_id >= 500 && parent_id < 600) || (parent_id >= 5000 && parent_id < 6000)) // b-hadrons
    return 2;
  if ((parent_id >= 400 && parent_id < 500) || (parent_id >= 4000 && parent_id < 5000)) // c-hadrons
    return 3;
  if ((parent_id >= 300 && parent_id < 400) || (parent_id >= 3000 && parent_id < 4000)) // s-hadrons
    return 4;

  return 5; // others
}

void MCNTupleFirstPass::FillSingleMuonParents(bool isMuon1){
// void FillSingleMuonParents(bool isMuon1, std::vector<int>& parent_barcodes, std::vector<int>& parent_ids){
  std::vector<int>::iterator itbar;
  std::vector<int>::iterator itprtbar;

  int ind = (isMuon1)? mpair->m1.ind : mpair->m2.ind;
  // std::cout << ind << ", " << mpair->m1.ind << ", " << mpair->m2.ind << std::endl;

  itbar = std::find(truth_barcode->begin(),truth_barcode->end(),ind);
  // std::cout << itbar - truth_barcode->begin() << ", " << truth_barcode->size() << std::endl;
  if (itbar != truth_barcode->end()){ // found the barcode of muon1 among all truth particles
    std::vector<int> parent_bars = truth_parents->at(itbar - truth_barcode->begin());
    if (isMuon1) mpair->m1_parent_barcode = parent_bars[0];
    else         mpair->m2_parent_barcode = parent_bars[0];

    // fill in histograms on number of parents & parent ids
    bool output_parents = false;

    if(parent_bars.size() > 1){
      h_numParents->Fill(parent_bars.size() - .5);
      output_parents = true;
    }
    else h_numParents->Fill(0.5);

    int parent_id;
    for (int parent_bar : parent_bars){
      itprtbar = std::find(truth_barcode->begin(),truth_barcode->end(),parent_bar);
      if (itprtbar != truth_barcode->end()){
        parent_id = abs(truth_id->at(itprtbar - truth_barcode->begin())) % 10000;
      }
      else{
        std::cout << "Error:: Barcode of the current parent not found among all truth particles, quitting" << std::endl;
        throw std::exception();
      }
    }

    if (isMuon1) mpair->m1_parent_id = parent_id;
    else         mpair->m2_parent_id = parent_id;

    if (parent_ids.size() == 1){
      int prt = ParentGrouping(parent_ids[0]);
      output_parents = (prt == -1);
      h_1parentGroups->Fill(prt - 0.5);
    }

    if (output_parents){
      for (int parent_id : parent_ids) m_unspecified_parent_file << parent_id << " ";
      m_unspecified_parent_file << std::endl;
    }

    ///////////////// first do single muons, get info, then 
    ///////////////// determine how we want to tag the single muons and group them
    ///////////////// Then we can do muon pairs /////////////////

    
  }else{
    std::cout << "Error:: Barcode of the current muon not found among all truth particles, quitting" << std::endl;
    throw std::exception();
  }

}

void MCNTupleFirstPass::FillMuonPairParents(){
  FillSingleMuonParents(true);
  FillSingleMuonParents(false);
  
  // change it up as a next step

}

void MCNTupleFirstPass::FillMuonPairTree(){
  int nsign = (mpair->same_sign)? 0:1;
  muonPairOutTree[nsign]->Fill();
  for (unsigned int idr = 0; idr < pms.ndRselcs; idr++){
    if (mpair->dr < pms.deltaR_thrsh[idr]){
      muonPairOutTreeBinned[idr][nsign]->Fill();
    }
  }
}


void MCNTupleFirstPass::ProcessData(){

  mpair = new MuonPairMC();
  // Long64_t nentries = 100000;
  Long64_t nentries = fChain->GetEntries();//number of events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {//loop over the events
  // for (Long64_t jentry=0; jentry<5;jentry++) {//loop over the events

    if(jentry%100000==0) cout<<"Processing "<<jentry<<" event out of "<<nentries<<" events"<<std::endl;

    int num_bytes = fChain->GetEntry(jentry);//read in an event
    if(num_bytes==0){
      std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
      throw std::exception();
    }
 

    //trigger requirement for event
    // if(!b_HLT_2mu4) continue;

    // std::vector<int> muon_index_list = {};
    // std::vector<int>::iterator it;

    int NPairs=muon_pair_muon1_pt->size();//number of muon pairs in the event
    h_numMuonPairsRaw->Fill(NPairs - 0.5);
    int NPairsAfter = 0;

    for(int i=0;i<NPairs;i++){//loop over all muon-pairs in the event

      // use ind to record barcode instead
      mpair->m1.ind     = static_cast<int>(muon_pair_muon1_bar->at(i));
      mpair->m2.ind     = static_cast<int>(muon_pair_muon2_bar->at(i));

      mpair->m1.pt    = fabs(muon_pair_muon1_pt->at(i))/1000.0;//pt of the first muon in the pair
      mpair->m2.pt    = fabs(muon_pair_muon2_pt->at(i))/1000.0;//pt of the second muon in the pair
      mpair->m1.eta   = muon_pair_muon1_eta->at(i);
      mpair->m2.eta   = muon_pair_muon2_eta->at(i);
      mpair->m1.phi   = muon_pair_muon1_phi->at(i);
      mpair->m2.phi   = muon_pair_muon2_phi->at(i);
      mpair->m1.charge  = static_cast<int>(muon_pair_muon1_ch->at(i));//sign of pt stores charge
      mpair->m2.charge  = static_cast<int>(muon_pair_muon2_ch->at(i));//sign of pt stores charge
      // mpair->m1.quality = muon_pair_muon1_quality->at(i);
      // mpair->m2.quality = muon_pair_muon2_quality->at(i);
      // mpair->m1.dP_overP = muon_deltaP_overP->at(mpair->m1.ind);
      // mpair->m2.dP_overP = muon_deltaP_overP->at(mpair->m2.ind);
      mpair->m1.ev_num = jentry;
      mpair->m2.ev_num = jentry;

      mpair->weight     = fabs(EventWeights->at(0)) * filter_effcy / nentries;
      mpair->minv       = truth_mupair_m->at(i)/1000.;
      mpair->pair_pt    = truth_mupair_pt->at(i)/1000;
      mpair->pair_y     = truth_mupair_y->at(i);
      mpair->asym       = truth_mupair_asym->at(i);
      mpair->acop       = truth_mupair_acop->at(i);

      // if(jentry%10000==0) std::cout << mpair->same_sign << std::endl;

      // if charge unequal gives 1 (opposite sign); otherwise gives 0 (same sign)
      h_cutAcceptance[mpair->m1.charge != mpair->m2.charge]->Fill(nocut + 0.5);
      //------------------------------------------------------------

      //Apply cuts

      //Trigger match for muon pair
      // if(dimuon_b_HLT_2mu4->at(i)==false) continue;

      if (!PassCuts())continue;
    
      //------------------------------------------------------------

      //Two things at this step: 
      //1) sort pt, eta, phi by pt
      //2) update the muon-pair values
      // mpair->Update();
      mpair->UpdateShort();

      // resonance cut
      if (IsResonance()) continue;

      // photo-production cut
      if (IsPhotoProduction()) continue;
      
      //------------------------------------------------------------

      // fill in muon parents (now that we have finished applying all the cuts)
      FillMuonPairParents();

      //------------------------------------------------------------

      NPairsAfter++;
      FillMuonPairTree(); // mode = 2 for MC truth

      // if(mode == 1){ // single muon information; no dR selection
      //   it = std::find(muon_index_list.begin(),muon_index_list.end(),mpair->m1.ind);
      //   if(it == muon_index_list.end()){ //muon1 index NOT found
      //     muon_index_list.push_back(mpair->m1.ind);
      //     tempmuon = &(mpair->m1);
      //     FillSingleMuonTree();
      //   }
      //   it = std::find(muon_index_list.begin(),muon_index_list.end(),mpair->m2.ind);
      //   if(it == muon_index_list.end()){ //muon1 index NOT found
      //     muon_index_list.push_back(mpair->m2.ind);
      //     tempmuon = &(mpair->m2);
      //     FillSingleMuonTree();
      //   }
      // }
      // else FillMuonPairTree(); // mode = 2
    }

    h_numMuonPairsAfter->Fill(NPairsAfter - 0.5);
  }//loop over events

  delete mpair;
}


void MCNTupleFirstPass::Run(){
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  if (mc_mode == "mc_truth_bb") filter_effcy = filter_effcy_bb;
  else if (mc_mode == "mc_truth_cc") filter_effcy = filter_effcy_cc;
  else{
      filter_effcy = -1000;
      std::cout << "MC mode has to be 'mc_truth_bb' or 'mc_truth_cc.' Else program will terminate." << std::endl;
  }
  std::cout << "The MC mode is " << mc_mode << ". Filter efficiency is " << filter_effcy << std::endl;

  // std::cout << "Mode = " << mode << ". Output file is " << m_outfile << std::endl;
  InitInput();
  InitOutput();
  ProcessData();
  m_unspecified_parent_file.close();
  // for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
  //   std::cout << "sign" << ksign << ", #entries: " << muonPairOutTree[ksign]->GetEntries() << std::endl;
  //   for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
  //     std::cout << "sign" << ksign << ", dr" << idr << ", #entries: " << muonPairOutTreeBinned[idr][ksign]->GetEntries() << std::endl;
  //   }
  // }

  // for loop over ibin: h->GetXaxis()->SetAxisLabel(ibin,label[ibin]);
  for (int ibin = 0; ibin < numCuts; ibin++){
    h_cutAcceptance[0]->GetXaxis()->SetBinLabel(ibin+1,cutLabels[ibin].c_str());
    h_cutAcceptance[1]->GetXaxis()->SetBinLabel(ibin+1,cutLabels[ibin].c_str());
  }

  for (int ibin = 0; ibin < nParentGroups; ibin++){
    h_1parentGroups->GetXaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
    h_2parentGroups->GetXaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
    h_2parentGroups->GetYaxis()->SetBinLabel(ibin+1,parentGroupLabels[ibin].c_str());
  }

  h_cutAcceptance[0]->Scale(1./h_cutAcceptance[0]->GetBinContent(1));
  h_cutAcceptance[1]->Scale(1./h_cutAcceptance[1]->GetBinContent(1));

  m_outfile->Write();

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "CPU time used (in seconds) for mode " << mode << " is " << cpu_time_used << std::endl;

}

#endif
