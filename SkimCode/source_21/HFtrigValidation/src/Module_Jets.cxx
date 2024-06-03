#include "HFtrigValidation/Module.h"
#include "HFtrigValidation/Module_Jets.h"
#include "xAODJet/JetContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include <TTree.h>





void Jets::Init(TTree *l_OutTree,int level_of_detail){
  m_level_of_detail=level_of_detail;

  char name[600];
  sprintf(name,"%s_pt"  ,m_BranchName.c_str());l_OutTree->Branch(name,&m_pt );
  sprintf(name,"%s_eta" ,m_BranchName.c_str());l_OutTree->Branch(name,&m_eta);
  sprintf(name,"%s_phi" ,m_BranchName.c_str());l_OutTree->Branch(name,&m_phi);
  sprintf(name,"%s_m"   ,m_BranchName.c_str());l_OutTree->Branch(name,&m_m  );
}

void Jets::InitTrackLinks(TTree *l_OutTree,std::string l_TrackParticleContainerKey){
  char name[600];
  sprintf(name,"%s_trk_index" ,m_BranchName.c_str());l_OutTree->Branch(name,&m_trk_index);
  m_TrackParticleContainerKey=l_TrackParticleContainerKey;
}

StatusCode Jets::ProcessTrackLinks(const xAOD::Jet *jet){
   static std::map<const xAOD::TrackParticle*,int> l_trk_index_temp;

   //first call in event is with jet=nullptr
   if(jet==nullptr){
     m_trk_index     .clear();
     l_trk_index_temp.clear();

     //retreieve the track container
     const xAOD::TrackParticleContainer *l_TrackParticleContainer;
     if(evtStore()->retrieve(l_TrackParticleContainer,m_TrackParticleContainerKey).isFailure()){
       //ATH_MSG_ERROR("Could not retrieve TrackParticleContainer with key "<<m_TrackParticleContainerKey.c_str());
       return StatusCode::FAILURE;
     }

     //index each track
     int trk_index=0;
     for(auto track:*l_TrackParticleContainer){
       l_trk_index_temp[track]=trk_index;
       trk_index++;
     }
   }
   //subsequent calls are with a jet sent in
   else{
     std::vector<int> indices_vector;

     const std::vector<ElementLink<DataVector<xAOD::IParticle>>> ptrackContainer=
     (jet->auxdata<std::vector<ElementLink<DataVector<xAOD::IParticle>>>>("constituentLinks"));

     for(auto& elink:ptrackContainer){
       if(elink.isValid()){
         const xAOD::TrackParticle *associated_track=dynamic_cast<const xAOD::TrackParticle*>(*elink);
         if(!associated_track){
           //ATH_MSG_ERROR();
           std::cout<<"ERROR in line "<<__LINE__<<" in file "<<__FILE__<<" :: nullptr track associated with jet of type "<<m_ContainerKey<<std::endl;
           //return StatusCode::FAILURE;
					 continue;
         }
         indices_vector.push_back(l_trk_index_temp[associated_track]);
       }
     }
     m_trk_index.push_back(indices_vector);
   }

   return StatusCode::SUCCESS;
}


StatusCode Jets::Process(){
   if(!IsEnabled()) return StatusCode::SUCCESS;

   m_pt.clear();
   m_eta.clear();
   m_phi.clear();
   m_m.clear();

   const xAOD::JetContainer* l_jet_container = nullptr;
   if(evtStore()->retrieve(l_jet_container,m_ContainerKey).isFailure()){
     std::cout<<"ERROR "<<__LINE__<<" "<<__FILE__<<"  Could not retrieve JetContainer with key "<<m_ContainerKey<<std::endl;
     return StatusCode::FAILURE;
   }

   if(m_TrackParticleContainerKey != "") CHECK(ProcessTrackLinks(nullptr));


   for(const auto *jet:(*l_jet_container)){
     m_pt.push_back (jet->pt ());
     m_eta.push_back(jet->eta());
     m_phi.push_back(jet->phi());
     m_m.push_back  (jet->m  ());

     if(m_TrackParticleContainerKey != "") CHECK(ProcessTrackLinks(jet));
   }

   return StatusCode::SUCCESS;
}

