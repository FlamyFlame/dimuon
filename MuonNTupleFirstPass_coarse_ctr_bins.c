#ifndef __MuonNTupleFirstPass_coarse_ctr_bins_C__
#define __MuonNTupleFirstPass_coarse_ctr_bins_C__

#include "MuonNTupleFirstPass_coarse_ctr_bins.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"

bool MuonNTupleFirstPass::PassCuts(){
  //Apply ALL CUTS but for resonances

  
  //require some quality cuts on the muons
  if((mpair->m1.quality&mpair->m2.quality&1  )==0) return false;//compair->m2ined muon
  if((mpair->m1.quality&mpair->m2.quality&8  )==0) return false;//Medium muon
  if((mpair->m1.quality&mpair->m2.quality&32 )==0) return false;//IDCuts
  if((mpair->m1.quality&mpair->m2.quality&256)==0) return false;//MuonCuts

  //cut on z0pair - order doesn't mpair->m1tter
  // double z0sinTheta1 = mpair->m1.z0*sin(2.0*atan(exp(-mpair->m1.eta)));
  // double z0sinTheta2 = mpair->m2.z0*sin(2.0*atan(exp(-mpair->m2.eta)));
  // float  z0_combined =sqrt(z0sinTheta1*z0sinTheta1 + z0sinTheta2*z0sinTheta2);
  // if(z0_combined>1.0) return false;

  if (fabs(mpair->m1.eta) > 2.4 || fabs(mpair->m2.eta) > 2.4) return false;
  if (mpair->m1.pt < 4 || mpair->m2.pt < 4) return false;
  // if( deltaP_overP_cut && (fabs(mpair->m1.dP_overP) > pms.deltaP_overP_thrsh || fabs(mpair->m2.dP_overP) > pms.deltaP_overP_thrsh) ) return false;
  if( fabs(mpair->m1.dP_overP) > pms.deltaP_overP_thrsh || fabs(mpair->m2.dP_overP) > pms.deltaP_overP_thrsh ) return false;
 
  return true;
}

bool MuonNTupleFirstPass::IsResonance(){
  // assuming we have already filled up the muon pair mpair
  if (!(mpair->same_sign)){ // opposite sign
    if (mpair->minv > pms.minv_upper) return true; // upper cut at 80 GeV

    bool isresonance = false;
    
    for (array<float,2> ires : pms.minv_cuts){
      if (mpair->minv > ires[0] && mpair->minv < ires[1]){
        isresonance = true;
        break;
      }
    }
      
    if (isresonance) return true;
  }

  return false; // same sign
}

bool MuonNTupleFirstPass::IsPhotoProduction(){
  return (!(mpair->same_sign) && mpair->asym < 0.05 && mpair->acop < 0.01);
}

void MuonNTupleFirstPass::FillSingleMuonTree(){
  muonOutTree->Fill();
  for (int jctr = 0; jctr < pms.nCtrBins; jctr++){
    if (tempmuon->ev_centrality >= 20*jctr && tempmuon->ev_centrality < 20*(jctr+1)){
      muonOutTreeBinned[jctr]->Fill();
    }
  }
}

void MuonNTupleFirstPass::FillMuonPairTree(){
  int nsign = (mpair->same_sign)? 0:1;
    muonPairOutTree[nsign]->Fill();
  for (unsigned int idr = 0; idr < pms.ndRselcs; idr++){
    if (mpair->dr < pms.deltaR_thrsh[idr]){
      for (unsigned int jctr = 0; jctr < pms.nCtrBins; jctr++){
        if (mpair->avg_centrality >= 20*jctr && mpair->avg_centrality < 20*(jctr+1)){
          muonPairOutTreeBinned[idr][jctr][nsign]->Fill();
        }
      }
    }
  }
}

void MuonNTupleFirstPass::ProcessData(){

  // Long64_t nentries = 100000;
  Long64_t nentries = fChain->GetEntries();//number of events
  for (Long64_t jentry=0; jentry<nentries;jentry++) {//loop over the events

    if(jentry%100000==0) cout<<"Processing "<<jentry<<" event out of "<<nentries<<" events"<<std::endl;

    int num_bytes = fChain->GetEntry(jentry);//read in an event
    if(num_bytes==0){
      std::cout<<"Error:: Read in event has size of zero bytes,  quitting"<<std::endl;
      throw std::exception();
    }
 

    //trigger requirement for event
    if(!b_HLT_mu4_mu4noL1) continue;


    std::vector<int> muon_index_list = {};
    std::vector<int>::iterator it;

    int NPairs=muon_pair_muon1_pt->size();//number of muon pairs in the event

    for(int i=0;i<NPairs;i++){//loop over all muon-pairs in the event

      mpair = new MuonPair();

      mpair->m1.ind     = muon_pair_muon1_index->at(i);
      mpair->m2.ind     = muon_pair_muon2_index->at(i);
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
      mpair->m1.ev_centrality = centrality;
      mpair->m2.ev_centrality = centrality;
      mpair->m1.ev_FCal_Et = FCal_Et;
      mpair->m2.ev_FCal_Et = FCal_Et;


      //------------------------------------------------------------

      //Apply cuts

      //Trigger match for muon pair
      if(dimuon_b_HLT_mu4_mu4noL1->at(i)==false) continue;

      if (!PassCuts())continue;
    
      //------------------------------------------------------------

      //Two things at this step: 
      //1) sort pt, eta, phi by pt
      //2) update the muon-pair values
      mpair->Update();

      // resonance cut
      if (IsResonance()) continue;

      // photo-production cut
      if (IsPhotoProduction()) continue;

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

void MuonNTupleFirstPass::Run(){
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  Init();
  InitOutput();
  ProcessData();
  m_outfile->Write();

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "CPU time used (in seconds) for mode " << mode << " is " << cpu_time_used << std::endl;

}

#endif
