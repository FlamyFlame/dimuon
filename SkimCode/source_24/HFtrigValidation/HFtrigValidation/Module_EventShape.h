#ifndef __MODULE_EVENTSHAPE_H__
#define __MODULE_EVENTSHAPE_H__

#include "HFtrigValidation/Module.h"
#include "Rtypes.h"

class EventShape:public Module{
  enum{
    NDET =2,//(FullCal, FCal)
      TOTAL_CALORIMETER  =0,
      FORWARD_CALORIMETER=1,
    NSIDE=3,//(+ve,-ve,combined)
      POSITIVE=0,
      NEGATIVE=1,
      COMBINED=2,
    NHAR=4,//V2-V3-v4-v5
  };

private:
  Float_t m_Calo_Et[NDET][NSIDE]      ;
  Float_t m_Calo_Qx[NDET][NSIDE][NHAR];
  Float_t m_Calo_Qy[NDET][NSIDE][NHAR];
  Float_t m_FCal_Et;
  Float_t m_FCal_Et_P;
  Float_t m_FCal_Et_N;
  UInt_t   m_year;
  Int_t   m_cent;

  int GetCentrality(float fcalET, const std::vector<float> & fCal_centr_boundaries);
  int GetCentralityPbPb2015(float FCal_Et);
  int GetCentralityPbPb2023(float FCal_Et);
  int GetCentralityPbPb2024(float FCal_Et);

public:
EventShape(ServiceHandle<StoreGateSvc> evtStore,
           std::string ContainerName,
           UInt_t Year=0,
           std::string BranchName="") 
        :m_year(Year), Module(evtStore,ContainerName,BranchName){};


  virtual void Init(TTree *l_OutTree,int level_of_detail=0) override;
  virtual StatusCode Process() override;
};
#endif
