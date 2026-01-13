#include "PowhegFullSimExtras.h"

template <class Derived>
void PowhegFullSimExtras<Derived>::InitInputExtra(){
    fChain->SetBranchAddress("muon_pair_muon1_pt"           , &muon_pair_muon1_pt);
    fChain->SetBranchAddress("muon_pair_muon1_eta"          , &muon_pair_muon1_eta);
    fChain->SetBranchAddress("muon_pair_muon1_phi"          , &muon_pair_muon1_phi);
    fChain->SetBranchAddress("muon_pair_muon1_ch"           , &muon_pair_muon1_ch);
    fChain->SetBranchAddress("muon_pair_muon1_bar"          , &muon_pair_muon1_bar);

    fChain->SetBranchAddress("muon_pair_muon2_pt"           , &muon_pair_muon2_pt);
    fChain->SetBranchAddress("muon_pair_muon2_eta"          , &muon_pair_muon2_eta);
    fChain->SetBranchAddress("muon_pair_muon2_phi"          , &muon_pair_muon2_phi);
    fChain->SetBranchAddress("muon_pair_muon2_ch"           , &muon_pair_muon2_ch);
    fChain->SetBranchAddress("muon_pair_muon2_bar"          , &muon_pair_muon2_bar);

    fChain->SetBranchStatus("muon_pair_muon1_pt"           ,1);
    fChain->SetBranchStatus("muon_pair_muon1_eta"              ,1);
    fChain->SetBranchStatus("muon_pair_muon1_phi"             ,1);
    fChain->SetBranchStatus("muon_pair_muon1_ch"             ,1);
    fChain->SetBranchStatus("muon_pair_muon1_bar"             ,1);

    fChain->SetBranchStatus("muon_pair_muon2_pt"             ,1);
    fChain->SetBranchStatus("muon_pair_muon2_eta"         ,1);
    fChain->SetBranchStatus("muon_pair_muon2_phi"              ,1);
    fChain->SetBranchStatus("muon_pair_muon2_ch"              ,1);
    fChain->SetBranchStatus("muon_pair_muon2_bar"              ,1);
}

template <class Derived>
bool PowhegFullSimExtras<Derived>::PassCutsExtra(){
}

template <class Derived>
void PowhegFullSimExtras<Derived>::FillMuonPairExtra(){
  	mpair->m1.ind     = muon_pair_muon1_bar->at(pair_ind);
  	mpair->m2.ind     = muon_pair_muon2_bar->at(pair_ind);
  	mpair->m1.charge  = muon_pair_muon1_ch->at(pair_ind);//sign of pt stores charge
  	mpair->m2.charge  = muon_pair_muon2_ch->at(pair_ind);//sign of pt stores charge
	
  	mpair->m1.pt    = fabs(muon_pair_muon1_pt->at(pair_ind))/1000.0;//pt of the first muon in the pair
  	mpair->m2.pt    = fabs(muon_pair_muon2_pt->at(pair_ind))/1000.0;//pt of the second muon in the pair
  	mpair->m1.eta   = muon_pair_muon1_eta->at(pair_ind);
  	mpair->m2.eta   = muon_pair_muon2_eta->at(pair_ind);
  	mpair->m1.phi   = muon_pair_muon1_phi->at(pair_ind);
  	mpair->m2.phi   = muon_pair_muon2_phi->at(pair_ind);
}