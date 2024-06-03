#ifndef __MODULE_JETS_H__
#define __MODULE_JETS_H__

#include "HFtrigValidation/Module.h"
#include "Rtypes.h"
#include "xAODJet/JetContainer.h"

class Jets:public Module{

private:
   std::vector<float> m_pt;
   std::vector<float> m_eta;
   std::vector<float> m_phi;
   std::vector<float> m_m;


   std::vector<std::vector<int>> m_trk_index;
   std::string m_TrackParticleContainerKey="";
   StatusCode ProcessTrackLinks(const xAOD::Jet* jet);

public:
Jets(ServiceHandle<StoreGateSvc> evtStore,
           std::string ContainerName,
           std::string BranchName="") 
        :Module(evtStore,ContainerName,BranchName){
  m_TrackParticleContainerKey="";
};


  virtual void Init(TTree *l_OutTree,int level_of_detail=0) override;
  virtual StatusCode Process() override;

   void InitTrackLinks(TTree *l_OutTree,std::string l_TrackParticleContainerKey);
};
#endif
