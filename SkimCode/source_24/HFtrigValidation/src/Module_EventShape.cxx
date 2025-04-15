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
  
  switch (m_year){
  case 2015:
    m_cent=GetCentralityPbPb2015(FCal_Et/1e6);
    break;
  case 2023:
    m_cent=GetCentralityPbPb2023(FCal_Et/1e6);
    break;
  case 2024:
    m_cent=GetCentralityPbPb2024(FCal_Et/1e6);
    break;
  default:
    m_cent=1000; // set to some nonsense value
  }

  return StatusCode::SUCCESS;
}

const std::vector<float> FCal_ET_Bins_PbPb2015={
4.26258 ,4.08137,3.91763 ,3.7635  ,3.61844 ,3.48077  ,3.34945 ,3.22397 ,3.10407 ,2.98931 ,//0-10
2.87864 ,2.77237 ,2.66999,2.57162 ,2.47658 ,2.38468  ,2.29572 ,2.21002 ,2.12711 ,2.04651 ,//10-20
1.96859 ,1.89316 ,1.81997 ,1.74932,1.68058 ,1.61434  ,1.55005 ,1.48744 ,1.42719 ,1.36875 ,//20-30
1.31197 ,1.25693 ,1.20373 ,1.15214 ,1.10211,1.05367  ,1.0068  ,0.961609,0.917795,0.87541 ,//30-40
0.834538,0.795018,0.756791,0.719896,0.684377,0.65018 ,0.617108,0.585275,0.554569,0.525092,//40-50
0.49675 ,0.46959 ,0.443549,0.418573,0.394518,0.371561,0.349697,0.328744,0.308686,0.289595,//50-60
0.27137 ,0.25407 ,0.237615,0.22199 ,0.207148,0.193096,0.179776,0.167193,0.155307,0.14414 ,//60-70
0.133573,0.123657,0.114352,0.105619,0.097388,0.089723,0.082548,0.075838,0.06956 ,0.063719,//70-80
0.0555  ,0.0505  ,0.0465  ,0.0425  ,0.0385  ,0.0355  ,0.0315  ,0.0275  ,0.0225   }; //80-89


const std::vector<float> FCal_ET_Bins_PbPb2023 = {
  4.51272, 4.32043, 4.15372, 3.99602, 3.84498, 3.69944, 3.55802, 3.42045, 3.28744, 3.15972, // 0-10
  3.03748, 2.92012, 2.80723, 2.69878, 2.59464, 2.49406, 2.39646, 2.3018, 2.21028, 2.12188, // 10-20
  2.03659, 1.95428, 1.87489, 1.79842, 1.72484, 1.65387, 1.58516, 1.51853, 1.45406, 1.39178, // 20-30
  1.33168, 1.27363, 1.21752, 1.16336, 1.11112, 1.06069, 1.0121, 0.965176, 0.919908, 0.876324, // 30-40
  0.834306, 0.793842, 0.754927, 0.71753, 0.681616, 0.647138, 0.61393, 0.581945, 0.55126, 0.52171, // 40-50
  0.493477, 0.466404, 0.440432, 0.415694, 0.392004, 0.369334, 0.347675, 0.327011, 0.30731, 0.288534, // 50-60
  0.270635, 0.253565, 0.237264, 0.221782, 0.207197, 0.19325, 0.179834, 0.167476, 0.155568, 0.144347, // 60-70
  0.13388, 0.12389, 0.114683, 0.105976, 0.0976472, 0.0901451, 0.082643, 0.0758922, 0.0695501, 0.063208, // 70-80
  0.0575959, 0.052731, 0.0478661, 0.0430012, 0.0388482 // 80-85
};

const std::vector<float> FCal_ET_Bins_PbPb2024 = {
  4.52146, 4.31811, 4.13353, 3.96046, 3.79717, 3.64263, 3.49592, 3.35557, 3.22165, 3.09343, // 0-10
  2.97016, 2.8517, 2.73778, 2.62777, 2.52149, 2.41863, 2.3191, 2.22271, 2.12924, 2.03863, // 10-20
  1.95086, 1.86565, 1.78306, 1.70296, 1.62528, 1.55017, 1.47738, 1.40701, 1.33893, 1.27333, // 20-30
  1.21019, 1.14936, 1.09078, 1.0345, 0.980452, 0.928742, 0.879123, 0.831534, 0.786079, 0.74255, // 30-40
  0.700943, 0.661196, 0.623259, 0.587086, 0.552536, 0.519523, 0.488175, 0.45821, 0.429779, 0.402687, // 40-50
  0.377083, 0.352811, 0.329731, 0.307802, 0.287017, 0.26733, 0.248698, 0.231103, 0.214482, 0.198731, // 50-60
  0.184059, 0.170295, 0.157174, 0.145029, 0.133661, 0.122843, 0.113043, 0.103605, 0.0950586, 0.0870635, // 60-70
  0.0792377, 0.0726958, 0.066154, 0.059703, 0.0546939, 0.0496847, 0.0446756, 0.0397848, 0.0365526, 0.0333205, // 70-80
  0.0300883, // 80-81
};

int EventShape::GetCentrality(float fcalET, const std::vector<float> & fCal_centr_boundaries) {
  // fCal_centr_boundaries is sorted descending
  // We want the first boundary that is <= fcalET.
  // That is exactly what lower_bound with std::greater does:
  auto it = std::lower_bound(fCal_centr_boundaries.begin(),
                             fCal_centr_boundaries.end(),
                             fcalET,
                             std::greater<float>());

  // centrality is the 0-based index
  int centrality = it - fCal_centr_boundaries.begin();

  if (centrality < 0) centrality = 0;
  if (centrality >= (int)fCal_centr_boundaries.size()) { // if centrality >= 85%, return -1
    return -1;
  }

  return centrality;
}


int EventShape::GetCentralityPbPb2015(float FCal_Et){
  return GetCentrality(FCal_Et, FCal_ET_Bins_PbPb2015);
}

int EventShape::GetCentralityPbPb2023(float FCal_Et){
  return GetCentrality(FCal_Et, FCal_ET_Bins_PbPb2023);
}

int EventShape::GetCentralityPbPb2024(float FCal_Et){
  return GetCentrality(FCal_Et, FCal_ET_Bins_PbPb2024);
}


