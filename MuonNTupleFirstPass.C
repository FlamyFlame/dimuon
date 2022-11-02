#ifndef __MuonNTupleFirstPass_C__
#define __MuonNTupleFirstPass_C__

#include "MuonNTupleFirstPass.h"
#include "Riostream.h"
#include "TLorentzVector.h"
#include "time.h"

bool MuonNTupleFirstPass::PassCuts(){
  //Apply ALL CUTS but for resonances

  bool pass = true;
  //require some quality cuts on the muons
  if((mpair->m1.quality&mpair->m2.quality&1  )==0) pass = false;//compair->m2ined muon
  if((mpair->m1.quality&mpair->m2.quality&8  )==0) pass = false;//Medium muon
  if((mpair->m1.quality&mpair->m2.quality&32 )==0) pass = false;//IDCuts
  if((mpair->m1.quality&mpair->m2.quality&256)==0) pass = false;//MuonCuts

  //cut on z0pair - order doesn't mpair->m1tter
  // double z0sinTheta1 = mpair->m1.z0*sin(2.0*atan(exp(-mpair->m1.eta)));
  // double z0sinTheta2 = mpair->m2.z0*sin(2.0*atan(exp(-mpair->m2.eta)));
  // float  z0_combined =sqrt(z0sinTheta1*z0sinTheta1 + z0sinTheta2*z0sinTheta2);
  // if(z0_combined>1.0) pass = false;

  if (fabs(mpair->m1.eta) > 2.4 || fabs(mpair->m2.eta) > 2.4) pass = false;
  if (mpair->m1.pt < 4 || mpair->m2.pt < 4) pass = false;
  // if( deltaP_overP_cut && (fabs(mpair->m1.dP_overP) > pms.deltaP_overP_thrsh || fabs(mpair->m2.dP_overP) > pms.deltaP_overP_thrsh) ) pass = false;
  if( fabs(mpair->m1.dP_overP) > pms.deltaP_overP_thrsh || fabs(mpair->m2.dP_overP) > pms.deltaP_overP_thrsh ) pass = false;
 
  return pass;
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

void MuonNTupleFirstPass::FillUnbinnedHistograms(){
  int nsign = (mpair->same_sign)? 0:1;

  for (unsigned int idr = 0; idr < pms.ndRselcs; idr++){
    if (mpair->dr < pms.deltaR_thrsh[idr]){
      h_FCal_Et[idr]->Fill(FCal_Et);
      h_pair_dP_overP[idr][nsign]->Fill(mpair->pair_dPoverP);
      h_Minv[idr][nsign]   ->Fill(mpair->minv);
      h_Dphi[idr][nsign]   ->Fill(mpair->dphi);
      h_Deta[idr][nsign]   ->Fill(mpair->deta);
      h_DR[idr][nsign]     ->Fill(mpair->dr);
      h_pt_lead[idr][nsign] ->Fill(mpair->m1.pt);
      h_eta_avg[idr][nsign] ->Fill(mpair->etaavg);
      h_pair_pt[idr][nsign]->Fill(mpair->pair_pt);
      h_pair_eta[idr][nsign]->Fill(mpair->pair_eta);
      h_pair_y[idr][nsign]->Fill(mpair->pair_y);
      h_eta_phi[idr][nsign]->Fill(mpair->phiavg,mpair->etaavg);
      h_eta1_eta2[idr][nsign]->Fill(mpair->m2.eta,mpair->m1.eta);
      h_pt1_pt2[idr][nsign]->Fill(mpair->m2.pt,mpair->m1.pt);
      h_eta_avg_Deta[idr][nsign]->Fill(mpair->deta,mpair->etaavg);        
      h_eta_avg_pair_eta[idr][nsign]->Fill(mpair->pair_eta,mpair->etaavg);
      h_ptlead_pair_pt[idr][nsign]->Fill(mpair->pair_pt,mpair->m1.pt);
    }
  }
}

void MuonNTupleFirstPass::FillCtrBinnedHistograms(){
  int nsign = (mpair->same_sign)? 0:1;

  for (unsigned int idr = 0; idr < pms.ndRselcs; idr++){
    if (mpair->dr < pms.deltaR_thrsh[idr]){
      for (unsigned int jctr = 0; jctr < pms.nCtrBins; jctr++){
        if (mpair->avg_centrality >= 20*jctr && mpair->avg_centrality < 20*(jctr+1)){
          h_ctrbin_pair_dP_overP[idr][jctr][nsign]->Fill(mpair->pair_dPoverP);
          h_ctrbin_Minv[idr][jctr][nsign]   ->Fill(mpair->minv);
          h_ctrbin_Dphi[idr][jctr][nsign]   ->Fill(mpair->dphi);
          h_ctrbin_Deta[idr][jctr][nsign]   ->Fill(mpair->deta);
          h_ctrbin_DR[idr][jctr][nsign]     ->Fill(mpair->dr);
          h_ctrbin_pt_lead[idr][jctr][nsign] ->Fill(mpair->m1.pt);
          h_ctrbin_eta_avg[idr][jctr][nsign] ->Fill(mpair->etaavg);
          h_ctrbin_pair_pt[idr][jctr][nsign]->Fill(mpair->pair_pt);
          h_ctrbin_pair_eta[idr][jctr][nsign]->Fill(mpair->pair_eta);
          h_ctrbin_pair_y[idr][jctr][nsign]->Fill(mpair->pair_y);
          h_ctrbin_eta_phi[idr][jctr][nsign]->Fill(mpair->phiavg,mpair->etaavg);
          h_ctrbin_eta1_eta2[idr][jctr][nsign]->Fill(mpair->m2.eta,mpair->m1.eta);
          h_ctrbin_pt1_pt2[idr][jctr][nsign]->Fill(mpair->m2.pt,mpair->m1.pt);
          h_ctrbin_eta_avg_Deta[idr][jctr][nsign]->Fill(mpair->deta,mpair->etaavg);        
          h_ctrbin_eta_avg_pair_eta[idr][jctr][nsign]->Fill(mpair->pair_eta,mpair->etaavg);
          h_ctrbin_ptlead_pair_pt[idr][jctr][nsign]->Fill(mpair->pair_pt,mpair->m1.pt);
        }
      }
    }
  }
}

void MuonNTupleFirstPass::FillPtBinnedHistograms(){

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
      if (mode != 1 && !(mpair->same_sign)){ //not raw (kill resonances) && opposite sign
        bool isresonance = false;
        for (array<float,2> ires : pms.minv_cuts){
          if (mpair->minv > ires[0] && mpair->minv < ires[1]) isresonance = true;
        }
        if (isresonance) continue;
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
      else{// only need pair information; does dR selection & binning

        if (mode == 1 || mode == 2){ // fill in unbinned hisograms
          //h_FCal_Et[idr]->Fill(FCal_Et);
          FillUnbinnedHistograms();
        }else if (mode == 4){
          FillMuonPairTree();
        }else if (mode == 5){
          FillCtrBinnedHistograms();
        }else if (mode == 6){
          FillPtBinnedHistograms();
        }else{
          std::cout << "The mode entered is undefined." << std::endl;
          std::cout << "Please enter a valid mode between 1 and 6." << std::endl;
          return;
        }
      }
    }
  }//loop over events
}

void MuonNTupleFirstPass::WriteOutput(){
  if (mode == 5){
    m_outfile->cd();
    // if (m_outfile->GetDirectory("ctr-binned") == nullptr)
    if (gDirectory->GetDirectory("ctr-binned") == nullptr){
      gDirectory->mkdir("ctr-binned");
    }
    gDirectory->cd("ctr-binned");
    gDirectory->Delete("h_*");

    for (unsigned int idr = 0; idr < pms.ndRselcs; idr++){
      for (unsigned int jctr = 0; jctr < pms.nCtrBins; jctr++){
        for (unsigned int ksign = 0; ksign < pms.nSigns; ksign++){
          h_ctrbin_pair_dP_overP[idr][jctr][ksign]->Write();
          h_ctrbin_Dphi[idr][jctr][ksign]->Write();
          h_ctrbin_Deta[idr][jctr][ksign]->Write();
          h_ctrbin_DR[idr][jctr][ksign]->Write();
          h_ctrbin_Minv[idr][jctr][ksign]->Write();
          h_ctrbin_pt_lead[idr][jctr][ksign]->Write();
          h_ctrbin_eta_avg[idr][jctr][ksign]->Write();
          h_ctrbin_pair_pt[idr][jctr][ksign]->Write();
          h_ctrbin_pair_eta[idr][jctr][ksign]->Write();
          h_ctrbin_pair_y[idr][jctr][ksign]->Write();
          h_ctrbin_eta_phi[idr][jctr][ksign]->Write();
          h_ctrbin_eta1_eta2[idr][jctr][ksign]->Write();
          h_ctrbin_eta_avg_Deta[idr][jctr][ksign]->Write();
          h_ctrbin_eta_avg_pair_eta[idr][jctr][ksign]->Write();
          h_ctrbin_pt1_pt2[idr][jctr][ksign]->Write();
          h_ctrbin_ptlead_pair_pt[idr][jctr][ksign]->Write();
        }
      }
    }
  }else if (mode == 6){
    m_outfile->cd();
    // if (m_outfile->GetDirectory("pt-binned") == nullptr){}
    if (gDirectory->GetDirectory("pt-binned") == nullptr){
      gDirectory->mkdir("pt-binned");
    }
    gDirectory->cd("pt-binned");
    gDirectory->Delete("h_*");

    for (unsigned int idr = 0; idr < pms.ndRselcs; idr++){
      for (unsigned int jpt = 0; jpt < pms.nPtBins; jpt++){
        for (unsigned int ksign = 0; ksign < pms.nSigns; ksign++){
          h_ptbin_pair_dP_overP[idr][jpt][ksign]->Write();
          h_ptbin_Dphi[idr][jpt][ksign]->Write();
          h_ptbin_Deta[idr][jpt][ksign]->Write();
          h_ptbin_DR[idr][jpt][ksign]->Write();
          h_ptbin_Minv[idr][jpt][ksign]->Write();
          h_ptbin_pt_lead[idr][jpt][ksign]->Write();
          h_ptbin_eta_avg[idr][jpt][ksign]->Write();
          h_ptbin_pair_pt[idr][jpt][ksign]->Write();
          h_ptbin_pair_eta[idr][jpt][ksign]->Write();
          h_ptbin_pair_y[idr][jpt][ksign]->Write();
          h_ptbin_eta_phi[idr][jpt][ksign]->Write();
          h_ptbin_eta1_eta2[idr][jpt][ksign]->Write();
          h_ptbin_eta_avg_Deta[idr][jpt][ksign]->Write();
          h_ptbin_eta_avg_pair_eta[idr][jpt][ksign]->Write();
          h_ptbin_pt1_pt2[idr][jpt][ksign]->Write();
          h_ptbin_ptlead_pair_pt[idr][jpt][ksign]->Write();
        }
      }
    }
  }
  else if(mode == 2){
    // for (unsigned int idr = 0; idr < pms.ndRselcs; idr++){
    //   std::cout << "h_FCal_Et[" << idr << "]\t" << h_FCal_Et[idr] << std::endl;
    //   for (unsigned int ksign = 0; ksign < pms.nSigns; ksign++){
    //     std::cout << "h_pair_dP_overP[" << idr << "]\t" << h_pair_dP_overP[idr][ksign] << std::endl;
    //     std::cout << "h_Dphi[" << idr << "]\t" << h_Dphi[idr][ksign] << std::endl;
    //     std::cout << "h_Deta[" << idr << "]\t" << h_Deta[idr][ksign] << std::endl;
    //     std::cout << "h_DR[" << idr << "]\t" << h_DR[idr][ksign] << std::endl;
    //     std::cout << "h_Minv[" << idr << "]\t" << h_Minv[idr][ksign] << std::endl;
    //     std::cout << "h_pt_lead[" << idr << "]\t" << h_pt_lead[idr][ksign] << std::endl;
    //     std::cout << "h_eta_avg[" << idr << "]\t" << h_eta_avg[idr][ksign] << std::endl;
    //     std::cout << "h_pair_pt[" << idr << "]\t" << h_pair_pt[idr][ksign] << std::endl;
    //     std::cout << "h_pair_eta[" << idr << "]\t" << h_pair_eta[idr][ksign] << std::endl;
    //     std::cout << "h_pair_y[" << idr << "]\t" << h_pair_y[idr][ksign] << std::endl;
    //     std::cout << "h_eta_phi[" << idr << "]\t" << h_eta_phi[idr][ksign] << std::endl;
    //     std::cout << "h_eta1_eta2[" << idr << "]\t" << h_eta1_eta2[idr][ksign] << std::endl;
    //     std::cout << "h_eta_avg_Deta[" << idr << "]\t" << h_eta_avg_Deta[idr][ksign] << std::endl;
    //     std::cout << "h_eta_avg_pair_eta[" << idr << "]\t" << h_eta_avg_pair_eta[idr][ksign] << std::endl;
    //     std::cout << "h_pt1_pt2[" << idr << "]\t" << h_pt1_pt2[idr][ksign] << std::endl;
    //     std::cout << "h_ptlead_pair_pt[" << idr << "]\t" << h_ptlead_pair_pt[idr][ksign] << std::endl;
    //   }
    // }

    m_outfile->cd();
    if (gDirectory->Get("h_*") != nullptr) gDirectory->Delete("h_*");

    for (unsigned int idr = 0; idr < pms.ndRselcs; idr++){
      h_FCal_Et[idr]->Write();
      // h_deltaP_overP[idr]->Write();
      
      for (unsigned int ksign = 0; ksign < pms.nSigns; ksign++){
        h_pair_dP_overP[idr][ksign]->Write();
        h_Dphi[idr][ksign]->Write();
        h_Deta[idr][ksign]->Write();
        h_DR[idr][ksign]->Write();
        h_Minv[idr][ksign]->Write();
        h_pt_lead[idr][ksign]->Write();
        h_eta_avg[idr][ksign]->Write();
        h_pair_pt[idr][ksign]->Write();
        h_pair_eta[idr][ksign]->Write();
        h_pair_y[idr][ksign]->Write();
        h_eta_phi[idr][ksign]->Write();
        h_eta1_eta2[idr][ksign]->Write();
        h_eta_avg_Deta[idr][ksign]->Write();
        h_eta_avg_pair_eta[idr][ksign]->Write();
        h_pt1_pt2[idr][ksign]->Write();
        h_ptlead_pair_pt[idr][ksign]->Write();
      }
    }
  }else{
    m_outfile->Write();
  }
}

void MuonNTupleFirstPass::Run(){
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  Init();
  InitOutput();
  InitHists();
  ProcessData();
  WriteOutput();

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "CPU time used (in seconds) for mode " << mode << " is " << cpu_time_used << std::endl;

}

#endif
