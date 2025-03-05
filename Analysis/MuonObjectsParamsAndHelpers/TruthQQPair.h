
#ifndef TruthQQPair_h
#define TruthQQPair_h

#include "Muon.h"
#include <iostream>
#include <string.h>
#include "TLorentzVector.h"
#include <math.h> 
#include <assert.h> 
#include <vector>


struct HQ{
    float pt;
    float eta;
    float phi;
};

class TruthQQPair{ 
// requirement: truth-level-MC QQbar pair with the same ancestors
//from hard scattering or single gluon splitting
public:
	HQ q1;
	HQ q2;
    int ev_num;
    float weight;
    int quark_type;
    float pt_lead;
    float pair_pt;
    float pair_y;
    float dpt;
    float deta;
    float dphi;
    float dr;
    float minv; // Q^2 for s channel production
    float asym; // asymmetry := (pT_lead - pT_sublead) / (pT_lead + pT_sublead)
    // float acop; // acoplanarity := (pi - |Dphi|) / pi
    int ancestor_group; // gg, gq, single g, qq
    // bool eta_same_sign;


    TruthQQPair(){quark_type = -1;}
    TruthQQPair(int quark);
    ~TruthQQPair(){}
    void Sort();
    void Update();
};


TruthQQPair::TruthQQPair(int quark){
    if (quark != 4 && quark != 5){
        std::cout << "Quark type needs to be 4 (c) or 5(b)! Quitting." << std::endl;
        throw std::exception();
    }
    quark_type = quark;
    q1.pt = 0;
    q2.pt = 0;
    q1.eta = 5;
    q2.eta = 5;
    q1.phi = 5;
    q2.phi = 5;
}



void TruthQQPair::Sort(){
    //sort pt, eta, phi by pt

    float temppt, tempeta, tempphi, tempd0, tempz0, tempdP_overP;
    int tempind, tempcharge, tempquality;
    if (q1.pt < q2.pt){
        temppt = q1.pt;
        tempphi = q1.phi;
        tempeta = q1.eta;
    
        q1.pt = q2.pt;
        q1.eta = q2.eta;
        q1.phi = q2.phi;
    
        q2.pt = temppt;
        q2.eta = tempeta;
        q2.phi = tempphi;
    }
}

void TruthQQPair::Update(){
    Sort();
    double PI=acos(-1.0);

    pt_lead = q1.pt;
    dpt = q1.pt - q2.pt;
    dphi = q1.phi - q2.phi;
    dphi = atan2(sin(dphi),cos(dphi));//fold dphi to [-pi,pi]
    deta = q1.eta - q2.eta;
    dr   = sqrt(dphi*dphi + deta*deta);
    // etaavg = (q1.eta + q2.eta)/2;
    // phiavg = (q1.phi + q2.phi)/2;

    TLorentzVector Q1, Q2, Q3;
    float quark_mass = (quark_type == 4)? 1.27 : 4.18;
    Q1.SetPtEtaPhiM(q1.pt,q1.eta,q1.phi,quark_mass);
    Q2.SetPtEtaPhiM(q2.pt,q2.eta,q2.phi,quark_mass);
    Q3=Q1+Q2;
    minv     = Q3.M();
    pair_pt  = Q3.Pt();
    pair_y   = Q3.Rapidity();

    asym = dpt / (q1.pt + q2.pt);
    assert (asym >= 0);
}

#endif
