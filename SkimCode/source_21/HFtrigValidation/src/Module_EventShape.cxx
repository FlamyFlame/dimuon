#include "HFtrigValidation/Module.h"
#include "HFtrigValidation/Module_EventShape.h"
#include "xAODHIEvent/HIEventShapeContainer.h"
#include <TTree.h>





void EventShape::Init(TTree *l_OutTree,int level_of_detail){
  m_level_of_detail=level_of_detail;

  char name[600];
  char branch_name[600];
  sprintf(name       ,"%sCalo_Et"              ,m_BranchName.c_str());
  sprintf(branch_name,"%sCalo_Et[%d][%d]/F"    ,m_BranchName.c_str(),NDET,NSIDE     );
  l_OutTree->Branch(name,&m_Calo_Et,branch_name);
  if(m_level_of_detail==0){
    sprintf(name       ,"%sCalo_Qx"              ,m_BranchName.c_str());
    sprintf(branch_name,"%sCalo_Qx[%d][%d][%d]/F",m_BranchName.c_str(),NDET,NSIDE,NHAR);
    l_OutTree->Branch(name,&m_Calo_Qx,branch_name);

    sprintf(name       ,"%sCalo_Qy"              ,m_BranchName.c_str());
    sprintf(branch_name,"%sCalo_Qy[%d][%d][%d]/F",m_BranchName.c_str(),NDET,NSIDE,NHAR);
    l_OutTree->Branch(name,&m_Calo_Qy,branch_name); 

  }
  else if(m_level_of_detail==1){
    sprintf(name       ,"%sCalo_Q"              ,m_BranchName.c_str());
    sprintf(branch_name,"%sCalo_Q[%d][%d][%d]/F",m_BranchName.c_str(),NDET,NSIDE,NHAR);
    l_OutTree->Branch(name,&m_Calo_Qx,branch_name);

    sprintf(name       ,"%sCalo_Psi"              ,m_BranchName.c_str());
    sprintf(branch_name,"%sCalo_Psi[%d][%d][%d]/F",m_BranchName.c_str(),NDET,NSIDE,NHAR);
    l_OutTree->Branch(name,&m_Calo_Qy,branch_name); 

  }
  sprintf(name       ,"%sFCal_Et"    ,m_BranchName.c_str());
  sprintf(branch_name,"%sFCal_Et/F"  ,m_BranchName.c_str());
  l_OutTree->Branch(name,&m_FCal_Et  ,branch_name); 

  sprintf(name       ,"%sFCal_Et_P"  ,m_BranchName.c_str());
  sprintf(branch_name,"%sFCal_Et_P/F",m_BranchName.c_str());
  l_OutTree->Branch(name,&m_FCal_Et_P,branch_name); 

  sprintf(name       ,"%sFCal_Et_N"  ,m_BranchName.c_str());
  sprintf(branch_name,"%sFCal_Et_N/F",m_BranchName.c_str());
  l_OutTree->Branch(name,&m_FCal_Et_N,branch_name); 

  sprintf(name       ,"%scentrality"  ,m_BranchName.c_str());
  sprintf(branch_name,"%scentrality/I",m_BranchName.c_str());
  l_OutTree->Branch(name,&m_cent,branch_name); 
}


StatusCode EventShape::Process(){
   if(!IsEnabled()) return StatusCode::SUCCESS;

   const xAOD::HIEventShapeContainer* l_HIevtShapeContainer;
   if(evtStore()->retrieve(l_HIevtShapeContainer,m_ContainerKey).isFailure()){
     std::cout<<" Could not retrieve HIEVentSHapeContainer with key "<<m_ContainerKey<<std::endl;
     return StatusCode::FAILURE;
   }


   double Calo_Et[NDET][NSIDE]      = {{0.0}};
   double Calo_Qx[NDET][NSIDE][NHAR]={{{0.0}}};
   double Calo_Qy[NDET][NSIDE][NHAR]={{{0.0}}};
   double FCal_Et                   =0.0;
   double FCal_Et_P                 =0.0;
   double FCal_Et_N                 =0.0;


   int size=l_HIevtShapeContainer->size();
   for(int ish=0;ish<size;ish++){
     const xAOD::HIEventShape *sh=l_HIevtShapeContainer->at(ish);
     float eta=(sh->etaMin()+sh->etaMax())*0.5;
     const std::vector<float> &c1=sh->etCos();
     const std::vector<float> &s1=sh->etSin();

     int layer=sh->layer();
     if(layer==21 || layer==22 || layer==23) {
       FCal_Et+=sh->et();
       if(eta>0) FCal_Et_P+=sh->et();
       if(eta<0) FCal_Et_N+=sh->et();
     }
     for(int idet=0;idet<NDET;idet++){
       for(int iside=0;iside<NSIDE;iside++){
         if(idet==FORWARD_CALORIMETER && !(layer==21 || layer==22 || layer==23)) continue;
         if(iside==POSITIVE   && (eta< 0) ) continue;
         if(iside==NEGATIVE   && (eta> 0) ) continue;
         Calo_Et[idet][iside] += sh->et();
         for(int ihar=0;ihar<NHAR;ihar++){
           Calo_Qx[idet][iside][ihar]+=c1[ihar+1];
           Calo_Qy[idet][iside][ihar]+=s1[ihar+1];
         }
       }
     }
   }


   for(int idet=0;idet<NDET;idet++){
     for(int iside=0;iside<NSIDE;iside++){
       for(int ihar=0;ihar<NHAR;ihar++){
         m_Calo_Et[idet][iside]      =Calo_Et[idet][iside];
         m_Calo_Qx[idet][iside][ihar]=Calo_Qx[idet][iside][ihar];
         m_Calo_Qy[idet][iside][ihar]=Calo_Qy[idet][iside][ihar];
         if(m_level_of_detail==1){
           double mag=sqrt(pow(Calo_Qx[idet][iside][ihar],2.0)+pow(Calo_Qy[idet][iside][ihar],2.0));
           double psi=atan2(Calo_Qy[idet][iside][ihar],Calo_Qx[idet][iside][ihar]);
           m_Calo_Qx[idet][iside][ihar]=mag;
           m_Calo_Qy[idet][iside][ihar]=psi;
         }
       }
     }
   }
   m_FCal_Et  =FCal_Et;
   m_FCal_Et_P=FCal_Et_P;
   m_FCal_Et_N=FCal_Et_N;
   m_cent=GetCentralityPbPb2015(FCal_Et/1e6);

   return StatusCode::SUCCESS;
}

float FCal_ET_Bins_DATA[100]={
4.26258 ,4.08137,3.91763 ,3.7635  ,3.61844 ,3.48077  ,3.34945 ,3.22397 ,3.10407 ,2.98931 ,//0-10
2.87864 ,2.77237 ,2.66999,2.57162 ,2.47658 ,2.38468  ,2.29572 ,2.21002 ,2.12711 ,2.04651 ,//10-20
1.96859 ,1.89316 ,1.81997 ,1.74932,1.68058 ,1.61434  ,1.55005 ,1.48744 ,1.42719 ,1.36875 ,//20-30
1.31197 ,1.25693 ,1.20373 ,1.15214 ,1.10211,1.05367  ,1.0068  ,0.961609,0.917795,0.87541 ,//30-40
0.834538,0.795018,0.756791,0.719896,0.684377,0.65018 ,0.617108,0.585275,0.554569,0.525092,//40-50
0.49675 ,0.46959 ,0.443549,0.418573,0.394518,0.371561,0.349697,0.328744,0.308686,0.289595,//50-60
0.27137 ,0.25407 ,0.237615,0.22199 ,0.207148,0.193096,0.179776,0.167193,0.155307,0.14414 ,//60-70
0.133573,0.123657,0.114352,0.105619,0.097388,0.089723,0.082548,0.075838,0.06956 ,0.063719,//70-80
0.0555  ,0.0505  ,0.0465  ,0.0425  ,0.0385  ,0.0355  ,0.0315  ,0.0275  ,0.0225  ,-1    }; //80-90 

int EventShape::GetCentralityPbPb2015(float FCal_Et){
  for(int i=0;i<90;i++){
    if(FCal_Et>FCal_ET_Bins_DATA[i]) return i;
  }
  return -1;
}


