#ifndef __MuonNTupleFirstPassPP_C__
#define __MuonNTupleFirstPassPP_C__

#include "MuonNTupleFirstPassPP.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"

bool MuonNTupleFirstPassPP::PassCuts(){
  //Apply ALL CUTS but for resonances

  
  //require some quality cuts on the muons
  if((mpair->m1.quality&mpair->m2.quality&1  )==0) return false;//compair->m2ined muon
  if (isTight){
    if ((mpair->m1.quality&mpair->m2.quality&16  )==0) return false;//Tight muon
  }else{
  if((mpair->m1.quality&mpair->m2.quality&8  )==0) return false;//Medium muon
  }
  if((mpair->m1.quality&mpair->m2.quality&32 )==0) return false;//IDCuts
  //if((mpair->m1.quality&mpair->m2.quality&256)==0) return false;//MuonCuts

  if (fabs(mpair->m1.eta) > 2.4 || fabs(mpair->m2.eta) > 2.4) return false;
  if (mpair->m1.pt < 4 || mpair->m2.pt < 4) return false;
  if( fabs(mpair->m1.dP_overP) > pms.deltaP_overP_thrsh || fabs(mpair->m2.dP_overP) > pms.deltaP_overP_thrsh ) return false;
 
  return true;
}

bool MuonNTupleFirstPassPP::IsResonance(){
  // assuming we have already filled up the muon pair mpair
  if (!(mpair->same_sign)){ // opposite sign
    if (mpair->minv > pms.minv_upper){ // upper cut at 60 GeV
      resonance_tagged_muon_index_list.push_back(mpair->m1.ind);
      resonance_tagged_muon_index_list.push_back(mpair->m2.ind);
      return true;
    }
    
    for (array<float,2> ires : pms.minv_cuts){
      if (mpair->minv > ires[0] && mpair->minv < ires[1]){
        resonance_tagged_muon_index_list.push_back(mpair->m1.ind);
        resonance_tagged_muon_index_list.push_back(mpair->m2.ind);
        return true;
      }
    }
  }

  return false; // same sign or opposite sign & not resonance
}

bool MuonNTupleFirstPassPP::IsPhotoProduction(){
  // A := |Delta PT| / (sum pT) < 0.05 && alpha := (pi-Dphi)/pi < 0.01
  // float A = (mpair->m1.pt - mpair->m2.pt) / (mpair->m1.pt + mpair->m2.pt);
  // assert (A >= 0);
  // float alpha = (pms.PI - fabs(mpair->dphi)) / pms.PI;
  return (!(mpair->same_sign) && mpair->asym < 0.05 && mpair->acop < 0.01);
}

void MuonNTupleFirstPassPP::FillSingleMuonTree(){
  muonOutTree->Fill();
}

void MuonNTupleFirstPassPP::FillMuonPairTree(){
  int nsign = (mpair->same_sign)? 0:1;
  // muonPairOutTree[nsign]->Fill();
  for (unsigned int idr = 0; idr < pms.ndRselcs; idr++){
    if (mpair->dr < pms.deltaR_thrsh[idr]){
      muonPairOutTreeBinned[idr][nsign]->Fill();
    }
  }
}


void MuonNTupleFirstPassPP::ProcessData(){

  // Long64_t nentries = 100000;
  Long64_t nentries = fChain->GetEntries();//number of events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {//loop over the events
  // for (Long64_t jentry=0; jentry<100000;jentry++) {//loop over the events

    if(jentry%100000==0) cout<<"Processing "<<jentry<<" event out of "<<nentries<<" events"<<std::endl;

    int num_bytes = fChain->GetEntry(jentry);//read in an event
    if(num_bytes==0){
      std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
      throw std::exception();
    }
 

    //trigger requirement for event
    if(!b_HLT_2mu4) continue;

    resonance_tagged_muon_index_list.clear(); // MUST CLEAR for each event!!
    std::vector<int> muon_index_list = {};
    std::vector<int>::iterator it;

    int NPairs=muon_pair_muon1_pt->size();//number of muon pairs in the event

    for(int i=0;i<NPairs;i++){//loop over all muon-pairs in the event

      mpair = new MuonPair();

      mpair->m1.ind     = muon_pair_muon1_index->at(i);
      mpair->m2.ind     = muon_pair_muon2_index->at(i);

      // // if m1 or m2 is resonance tagged, do not record the pair
      // std::string ss_str = (muon_pair_muon1_pt ->at(i) * muon_pair_muon2_pt ->at(i) > 0)? ", same" : ", opp";
      // it = std::find(resonance_tagged_muon_index_list.begin(),resonance_tagged_muon_index_list.end(),mpair->m1.ind);
      // if(it != resonance_tagged_muon_index_list.end()) {
      //   // cout << jentry << ", " << std::to_string(i) << "-th pair, " << mpair->m1.ind << ss_str << endl;
      //   continue;
      // }

      // it = std::find(resonance_tagged_muon_index_list.begin(),resonance_tagged_muon_index_list.end(),mpair->m2.ind);
      // if(it != resonance_tagged_muon_index_list.end()) {
      //   // cout << jentry << ", " << std::to_string(i) << "-th pair, " << mpair->m2.ind << ss_str << endl;
      //   continue;
      // }

      mpair->m1.pt    = fabs(muon_pair_muon1_pt->at(i))/1000.0;//pt of the first muon in the pair
      mpair->m2.pt    = fabs(muon_pair_muon2_pt->at(i))/1000.0;//pt of the second muon in the pair
      mpair->m1.eta   = muon_pair_muon1_eta->at(i);
      mpair->m2.eta   = muon_pair_muon2_eta->at(i);
      mpair->m1.phi   = muon_pair_muon1_phi->at(i);
      mpair->m2.phi   = muon_pair_muon2_phi->at(i);
      mpair->m1.d0    = muon_pair_muon1_d0 ->at(i);
      mpair->m2.d0    = muon_pair_muon2_d0 ->at(i);
      mpair->m1.z0    = muon_pair_muon1_z0 ->at(i);
      mpair->m2.z0    = muon_pair_muon2_z0 ->at(i);
      mpair->m1.charge  =(muon_pair_muon1_pt ->at(i) > 0)? 1:-1;//sign of pt stores charge
      mpair->m2.charge  =(muon_pair_muon2_pt ->at(i) > 0)? 1:-1;//sign of pt stores charge
      mpair->m1.quality = muon_pair_muon1_quality->at(i);
      mpair->m2.quality = muon_pair_muon2_quality->at(i);
      mpair->m1.dP_overP = muon_deltaP_overP->at(mpair->m1.ind);
      mpair->m2.dP_overP = muon_deltaP_overP->at(mpair->m2.ind);
      mpair->m1.ev_num = jentry;
      mpair->m2.ev_num = jentry;


      //------------------------------------------------------------

      //Apply cuts

      //Trigger match for muon pair
      if(dimuon_b_HLT_2mu4->at(i)==false) continue;

      if (!PassCuts())continue;
    
      //------------------------------------------------------------

      //Two things at this step: 
      //1) sort pt, eta, phi by pt
      //2) update the muon-pair values
      mpair->Update();

      // resonance cut
      if (IsResonance()) continue;

      // photo-production cut
      if (IsPhotoProduction()){
        // h_dphi_failing_photoprod->Fill(mpair->dphi);
        // h_asym_acop_failing_photoprod->Fill(mpair->acop,mpair->asym);
        continue;
      }
      
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


void MuonNTupleFirstPassPP::Run(){
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  std::cout << "Mode = " << mode << ". Output file is " << m_outfile << std::endl;
  InitInput();
  InitOutput();
  ProcessData();
  // for (unsigned int ksign = 0; ksign < ParamsSet::nSigns; ksign++){
  //   std::cout << "sign" << ksign << ", #entries: " << muonPairOutTree[ksign]->GetEntries() << std::endl;
  //   for (unsigned int idr = 0; idr < ParamsSet::ndRselcs; idr++){
  //     std::cout << "sign" << ksign << ", dr" << idr << ", #entries: " << muonPairOutTreeBinned[idr][ksign]->GetEntries() << std::endl;
  //   }
  // }
  m_outfile->Write();

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "CPU time used (in seconds) for mode " << mode << " is " << cpu_time_used << std::endl;

}

#endif
