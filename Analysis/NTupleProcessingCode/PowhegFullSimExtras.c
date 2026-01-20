#include "PowhegFullSimExtras.h"

template <class Derived>
void PowhegFullSimExtras<Derived>::InitInputExtra(){
    fChain()->SetBranchAddress("muon_pair_muon1_pt"           , &muon_pair_muon1_pt);
    fChain()->SetBranchAddress("muon_pair_muon1_eta"          , &muon_pair_muon1_eta);
    fChain()->SetBranchAddress("muon_pair_muon1_phi"          , &muon_pair_muon1_phi);
    fChain()->SetBranchAddress("muon_pair_muon1_index"          , &muon_pair_muon1_index);

    fChain()->SetBranchAddress("muon_pair_muon2_pt"           , &muon_pair_muon2_pt);
    fChain()->SetBranchAddress("muon_pair_muon2_eta"          , &muon_pair_muon2_eta);
    fChain()->SetBranchAddress("muon_pair_muon2_phi"          , &muon_pair_muon2_phi);
    fChain()->SetBranchAddress("muon_pair_muon2_index"          , &muon_pair_muon2_index);

    fChain()->SetBranchStatus("muon_pair_muon1_pt"           ,1);
    fChain()->SetBranchStatus("muon_pair_muon1_eta"              ,1);
    fChain()->SetBranchStatus("muon_pair_muon1_phi"             ,1);
    fChain()->SetBranchStatus("muon_pair_muon1_ch"             ,1);
    fChain()->SetBranchStatus("muon_pair_muon1_index"             ,1);

    fChain()->SetBranchStatus("muon_pair_muon2_pt"             ,1);
    fChain()->SetBranchStatus("muon_pair_muon2_eta"         ,1);
    fChain()->SetBranchStatus("muon_pair_muon2_phi"              ,1);
    fChain()->SetBranchStatus("muon_pair_muon2_ch"              ,1);
    fChain()->SetBranchStatus("muon_pair_muon2_index"              ,1);
}

template <class Derived>
void PowhegFullSimExtras<Derived>::FillMuonPairExtra(int pair_ind){
    mpair()->m1.ind     = muon_pair_muon1_index->at(pair_ind);
    mpair()->m2.ind     = muon_pair_muon2_index->at(pair_ind);
    
    mpair()->m1.pt    = fabs(muon_pair_muon1_pt->at(pair_ind))/1000.0;//pt of the first muon in the pair
    mpair()->m2.pt    = fabs(muon_pair_muon2_pt->at(pair_ind))/1000.0;//pt of the second muon in the pair
    mpair()->m1.eta   = muon_pair_muon1_eta->at(pair_ind);
    mpair()->m2.eta   = muon_pair_muon2_eta->at(pair_ind);
    mpair()->m1.phi   = muon_pair_muon1_phi->at(pair_ind);
    mpair()->m2.phi   = muon_pair_muon2_phi->at(pair_ind);

    mpair()->m1.d0    = muon_pair_muon1_d0 ->at(pair_ind);
    mpair()->m2.d0    = muon_pair_muon2_d0 ->at(pair_ind);
    mpair()->m1.z0    = muon_pair_muon1_z0 ->at(pair_ind);
    mpair()->m2.z0    = muon_pair_muon2_z0 ->at(pair_ind);
    mpair()->m1.charge  =(muon_pair_muon1_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge
    mpair()->m2.charge  =(muon_pair_muon2_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge
    mpair()->m1.quality = muon_pair_muon1_quality->at(pair_ind);
    mpair()->m2.quality = muon_pair_muon2_quality->at(pair_ind);
    mpair()->m1.dP_overP = muon_deltaP_overP->at(mpair()->m1.ind);
    mpair()->m2.dP_overP = muon_deltaP_overP->at(mpair()->m2.ind);
    
    mpair()->m1.trk_pt      = fabs(muon_pair_muon1_trk_pt->at(pair_ind))/1000.0;//pt of the first muon in the pair
    mpair()->m2.trk_pt      = fabs(muon_pair_muon2_trk_pt->at(pair_ind))/1000.0;//pt of the second muon in the pair
    mpair()->m1.trk_eta     = muon_pair_muon1_trk_eta->at(pair_ind);
    mpair()->m2.trk_eta     = muon_pair_muon2_trk_eta->at(pair_ind);
    mpair()->m1.trk_phi     = muon_pair_muon1_trk_phi->at(pair_ind);
    mpair()->m2.trk_phi     = muon_pair_muon2_trk_phi->at(pair_ind);

    if (turn_on_track_charge){
        mpair()->m1.trk_charge  =(muon_pair_muon1_trk_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge
        mpair()->m2.trk_charge  =(muon_pair_muon2_trk_pt ->at(pair_ind) > 0)? 1:-1;//sign of pt stores charge    
    }else{ // do not turn on track charge: set to nonsense
        mpair()->m1.trk_charge  = 0;
        mpair()->m2.trk_charge  = 0;
    }
}