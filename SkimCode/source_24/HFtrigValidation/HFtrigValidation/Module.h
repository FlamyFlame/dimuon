#ifndef __MODULE_H__
#define __MODULE_H__

#include "AthenaBaseComps/AthAlgorithm.h"
#include "string"


class TTree;

class Module{
  protected:
  ServiceHandle<StoreGateSvc> m_evtStore;
  ServiceHandle<StoreGateSvc>& evtStore(){return m_evtStore;}

  std::string m_ContainerKey;
  std::string m_BranchName;
  bool m_is_enabled=false;
  int  m_level_of_detail=0;

  bool IsEnabled(){return m_is_enabled;}

  public:
  Module(ServiceHandle<StoreGateSvc> evtStore,
         std::string ContainerName,
         std::string BranchName="")
         :m_evtStore(evtStore){
     m_ContainerKey=ContainerName;
     m_BranchName=BranchName;
     if(m_ContainerKey!="") m_is_enabled=true;
  };

  virtual void Init(TTree *l_OutTree,int level_of_detail=0)=0;
  virtual StatusCode Process()=0; 
};
#endif
