#include "HFtrigValidation/TrigRates.h"
#include "HFtrigValidation/MyUtils.h"
#include "HFtrigValidation/Module_EventShape.h"
#include "HFtrigValidation/Module_Jets.h"


#include <GaudiKernel/ITHistSvc.h>
#include <GaudiKernel/ServiceHandle.h>

#include "xAODHIEvent/HIEventShapeContainer.h"
#include "xAODTrigger/TrigDecision.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTrigger/EnergySumRoI.h"
#include "xAODForward/ZdcModuleContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertexContainer.h"

#include "GeneratorObjects/McEventCollection.h"

#include <TTree.h>
#include <cmath>
#include <Riostream.h>

enum Track{
  StoreBasic    =2,
  StoreDetails  =4,
  StoreTruthLink=8,
};

enum Truth{
  Basic        =1,
  StoreParents =2,
  StoreMuSingle=4,
  StoreMuPairs =8,
};

bool TrigRates::IsPrimaryParticle(const xAOD::TruthParticle* particle){
  if(!particle) return false;
  if(particle->status()!=1                             ) return false;
  if(particle->barcode()>=2e5 || particle->barcode()==0) return false;
  if(particle->charge()==0                             ) return false;
  return true;
}



TrigRates::TrigRates(const std::string& name, ISvcLocator* pSvcLocator)
  : AthAlgorithm(name,pSvcLocator),
  m_trigTool("Trig::TrigDecisionTool/TrigDecisionTool") 
{
  declareProperty("IsEvgen"               ,  m_is_evgen                    =false         );//true/false
  declareProperty("UseGRL"                ,  m_use_GRL                     =false         );//true/false
  declareProperty("MaxZvtx"               ,  m_MaxZvtx                     =250           );
  
  declareProperty("RunYear"               ,  m_year                        =0          );

  declareProperty("StoreAllEvents"        ,  m_StoreAllEvents              =true          );//true,false
  declareProperty("UseTrigger"            ,  m_use_trigger                 =true          );//true,false
  declareProperty("StoreL1Decision"       ,  m_store_L1                    =0             );//set to 1 for storing L1 TBP/TAP/TAV values 
  declareProperty("TriggerChains"         ,  m_Trigger_Chains              ="HLT_.*"      );
  declareProperty("MuonTriggerChains"     ,  m_Muon_Trigger_Chains         ="HLT_mu4"     );
  declareProperty("DiMuonTriggerChains"   ,  m_DiMuon_Trigger_Chains       ="HLT_2mu4"    );

  declareProperty("EventInfoKey"          ,  m_EventInfo_key               ="EventInfo"          );
  declareProperty("TrkContainerKey"       ,  m_trk_container_key           ="InDetTrackParticles");
  declareProperty("PixContainerKey"       ,  m_pix_container_key           ="InDetPixelTrackParticles");
  declareProperty("VtxContainerKey"       ,  m_vtx_container_key           ="PrimaryVertices"    );
  declareProperty("L1TEContainerKey"      ,  m_L1TE_container_key          ="LVL1EnergySumRoI"   );
  declareProperty("MuonsKey"              ,  m_muons_key                   ="Muons"              );
  declareProperty("HIEventShapeContainerKey",  m_HIEventShapeContainer_key ="HIEventShape"       );
  declareProperty("TruthVtxContainerKey"  ,  m_truth_vtx_container_key     ="TruthVertices"      );
  declareProperty("TruthContainerKey"     ,  m_truth_container_key         ="TruthParticles"     );
  declareProperty("TrackJetContainerKeys" ,  m_track_jet_container_keys    ={}                   );
  declareProperty("METContainerKey"       ,  m_met_container_key           ="MET_Calo"           );
  declareProperty("HLTMuonsKey"           ,  m_hlt_muons_key               ="HLT_MuonsCB_RoI"    );

  declareProperty("StoreEventInfo"        ,  m_store_EventInfo             =1             );//bitflag
  declareProperty("StoreTracks"           ,  m_store_tracks                =1             );//bitmap
  declareProperty("StorePixTracks"        ,  m_store_pix_tracks            =1             );//bitmap
  declareProperty("StoreVtx"              ,  m_store_Vtx                   =true          );//true/false
  declareProperty("StoreL1TE"             ,  m_store_L1TE                  =true          );//true/false

  declareProperty("StoreSingleMuon"       ,  m_store_single_muon           =true          );//true/false
  declareProperty("StoreAcoplanarMuon"    ,  m_store_acoplanar_muon        =true          );//true/false
  declareProperty("StoreMuonTruth"        ,  m_store_muon_truth            =true          );//true/false

  declareProperty("StoreTruthVtx"         ,  m_store_truth_Vtx             =false         );//true/false
  declareProperty("StoreTruth"            ,  m_store_truth                 =0             );//bitflag
  declareProperty("TruthMinPT"            ,  m_min_pT_Truth                =500           );//Min Truth Partilce pT in MeV
  declareProperty("ApplyMuonCalibrations" ,  m_ApplyMuonCalibrations       =false         );//Muon Calibrations
  declareProperty("StoreMcTruth"          ,  m_store_McEvent               =0             );

  declareProperty("TrackSelectionTool_HITight"    , m_trkSelTool_HITight                    );
  declareProperty("TrackSelectionTool_HILoose"    , m_trkSelTool_HILoose                    );
  declareProperty("TrackSelectionTool_MinBias"    , m_trkSelTool_MinBias                    );

  declareProperty("IsZdcCalib"            , m_is_Zdc_Calib                 =false          );//true only for Zdc Calib Stream
  declareProperty("StoreZdc"              , m_store_Zdc                    =0           );//bitflag
  declareProperty("ZdcAuxSuffix"          , m_ZdcAuxSuffix                 =""          );//AuxSuffix for Zdc Reprocessing

  declareProperty("IsRun3"                , m_isRun3 = true                 );

  declareProperty("MuonSelectionTool"     , m_muonSelection                 );
  #if defined(__ATHENA_21p2__) || defined(__ATHENA_24p2__)
    declareProperty("TriggerMatchTool"    , m_matchTool                     );
  #else
    declareProperty("MuonTriggerMatchTool", m_matchTool                     );
  #endif
  declareProperty("GRLTool"               , m_grlTool                       );
  declareProperty("muonCorrectionTool"    , m_muonCalibrationAndSmearingTool);
  declareProperty("use_effi_SF_tool"      , m_use_effi_SF_tool           );
  declareProperty("effi_SF_tool_medium"   , m_effi_SF_tool_medium           );
  declareProperty("effi_SF_tool_tight"    , m_effi_SF_tool_tight            );
  #ifndef __ATHENA_24p2__
    declareProperty("ZDCAnalysisTool"       , m_ZDCAnalysisTool               );
    declareProperty("HIPileupTool"          , m_HIPileupTool                  );
  #endif
}



StatusCode TrigRates::initialize(){
   ATH_MSG_INFO("Starting initialize()");


   //In this case we only setup the Truth branches
   if(m_is_evgen==true){
     ServiceHandle<ITHistSvc> histSvc("THistSvc",name());
     CHECK( histSvc.retrieve() );

     m_OutTree=new TTree("HeavyIonD3PD","HeavyIonD3PD");
     CHECK(histSvc->regTree("/MYSTREAM/HeavyIonD3PD",m_OutTree));

     if(m_store_truth_Vtx   ) InitTruthVertex(m_OutTree);
     if(m_store_truth       ) InitTruth      (m_OutTree);
     if(m_store_McEvent     ) InitMcEvents   (m_OutTree);
     return StatusCode::SUCCESS;
   }

//retreive the tools //THIS IS OPTIONAL!
//-------------------------------------------------------------------------
   if(m_use_GRL      ) CHECK( m_grlTool .retrieve() );

   if(m_use_trigger  ){
     if(m_isRun3){
       #ifdef __ATHENA_24p2__
       //https://twiki.cern.ch/twiki/bin/view/Atlas/R22TriggerAnalysis
       m_trigTool.setTypeAndName("Trig::TrigDecisionTool/TrigDecisionTool");
       CHECK(m_trigTool.setProperty("NavigationFormat","TrigComposite"));
       CHECK(m_trigTool.setProperty("HLTSummary","HLTNav_Summary_AODSlimmed"));
       CHECK( m_trigTool.initialize() );
       //CHECK( m_trigTool.retrieve() );
       #endif
     }
     else{
       //CHECK(m_trigTool.setProperty("NavigationFormat","TriggerElement"));
       CHECK( m_trigTool.retrieve() );
     }
   }


   if(m_muons_key!="" &&(m_store_single_muon || m_store_acoplanar_muon) && m_use_trigger ){
     CHECK( m_matchTool.retrieve());
     //CHECK( m_matchTool.setProperty("TrigDecisionTool",m_trigTool.getHandle()));
     if(m_ApplyMuonCalibrations){
       CHECK( m_muonCalibrationAndSmearingTool.retrieve());
       CHECK( m_muonCalibrationAndSmearingTool->initialize());
     }
     if (m_use_effi_SF_tool){
       CHECK(m_effi_SF_tool_medium .retrieve  ());
       CHECK(m_effi_SF_tool_tight  .retrieve  ());
       CHECK(m_effi_SF_tool_medium->initialize());
       CHECK(m_effi_SF_tool_tight ->initialize());      
     }
   }
   CHECK(m_trkSelTool_HITight.retrieve() );
   CHECK(m_trkSelTool_HILoose.retrieve() );
   CHECK(m_trkSelTool_MinBias.retrieve() );

   #ifndef __ATHENA_24p2__
     CHECK(m_ZDCAnalysisTool.retrieve   ());
     CHECK(m_ZDCAnalysisTool->initialize());

     CHECK(m_HIPileupTool   .retrieve   ());
     //std::cout<<"retrieve done"<<std::endl;
     //CHECK(m_HIPileupTool   ->initialize());
     //std::cout<<"initialize done"<<std::endl;
   #endif
//-------------------------------------------------------------------------



//Create and Book OutPutTree
//Trigger branches are added later in first loop of execute()
//-------------------------------------------------------------------------
   ServiceHandle<ITHistSvc> histSvc("THistSvc",name()); 
   CHECK( histSvc.retrieve() );

   m_OutTree=new TTree("HeavyIonD3PD","HeavyIonD3PD");
   CHECK(histSvc->regTree("/MYSTREAM/HeavyIonD3PD",m_OutTree));

   if(m_store_Zdc         ) InitZdc   (m_OutTree);
   if(m_store_EventInfo>0 ) InitEventInfo(m_OutTree);

   if(m_is_Zdc_Calib==false){
     if(m_store_tracks      ) InitTracks   (m_OutTree);
     if(m_store_pix_tracks  ) InitPixTracks(m_OutTree);
     if(m_store_Vtx         ) InitVertex   (m_OutTree);
     if(m_store_L1TE        ) InitL1TE     (m_OutTree);
     if(m_muons_key!="" &&(m_store_single_muon || m_store_acoplanar_muon)) InitMuons    (m_OutTree);
     if(m_store_truth_Vtx   ) InitTruthVertex(m_OutTree);
     if(m_store_truth       ) InitTruth      (m_OutTree);
     if(m_met_container_key!="") InitMET        (m_OutTree); 
  
  
     Module* evtshape_module;
     if(m_HIEventShapeContainer_key != ""){
       evtshape_module=new EventShape(evtStore(),m_HIEventShapeContainer_key, m_year);
       evtshape_module->Init(m_OutTree,1); 
       m_modules.push_back(evtshape_module);
     }
  
     for(auto jet_container:m_track_jet_container_keys){
       Module* Jet_Module=new Jets(evtStore(), jet_container, jet_container);
       Jet_Module->Init(m_OutTree);
       ((Jets*)Jet_Module)->InitTrackLinks(m_OutTree, m_trk_container_key);
       m_modules.push_back(Jet_Module);
     }
   }
//-------------------------------------------------------------------------


   return StatusCode::SUCCESS;
}







StatusCode TrigRates::execute(){
   ATH_MSG_DEBUG("Starting execute()");

   static int ievent=0;
   if(ievent%100==0) std::cout<<"Running Event "<<ievent<<std::endl;
   ievent++;

   if(m_use_trigger && m_is_evgen==false){
     static bool b_AddTriggerBranches=false;
     if(b_AddTriggerBranches==false) AddTriggerBranches();
     b_AddTriggerBranches=true;
   }



   bool pass_cleaning_cuts=true;
   if(m_is_evgen==false) {CHECK(CleaningCuts(pass_cleaning_cuts));}
   if(!pass_cleaning_cuts){
     ATH_MSG_DEBUG("Failed Cleaning cuts");
     return StatusCode::SUCCESS;
   }
   else ATH_MSG_DEBUG("Passed Cleaning cuts");


   if(m_use_trigger && m_is_evgen==false){
     bool pass_trigger_cuts =false;
     CHECK(ProcessTriggers(pass_trigger_cuts));
     //keep event even if it didnot pass any triggers if m_StoreAllEvents==true
     if(m_StoreAllEvents==false && pass_trigger_cuts==false){
       return StatusCode::SUCCESS;
     }
   }

   //Clear trigger matching
   for(auto trigger_match:m_trigger_match_map   ) trigger_match.second->clear();
   for(auto trigger_match:m_trigger_match_map_V2) trigger_match.second->clear();
   for(auto trigger_match:m_trigger_match_map_V3) trigger_match.second->clear();

   if(m_is_Zdc_Calib==true){
     if(m_store_EventInfo>0)    CHECK(ProcessEventInfo   ());
     if(m_store_Zdc        )    CHECK(ProcessZdc         ());
   }
   else{
     if(m_is_evgen==false){
       if(m_store_EventInfo>0)    CHECK(ProcessEventInfo   ());
       if(m_store_tracks     )    CHECK(ProcessTracks      ());
       if(m_store_pix_tracks )    CHECK(ProcessPixTracks   ());
       if(m_store_Vtx        )    CHECK(ProcessVertex      ());
       if(m_store_L1TE       )    CHECK(ProcessL1TE        ());
       //The lower line must be run after ProcessTracks for pointer to track contaainer to be valid
       if(m_muons_key!="" &&(m_store_single_muon || m_store_acoplanar_muon))  CHECK(ProcessMuons());
       for(auto module:m_modules) CHECK(module->Process()); 
       if(m_met_container_key!="")    CHECK(ProcessMET         ());
     }
     if(m_store_truth_Vtx  )    CHECK(ProcessTruthVertex ());
     if(m_store_truth      )    CHECK(ProcessTruth       ());//must be run after "ProcessTracks()"
     if(m_store_McEvent    )    CHECK(ProcessMcEvents    ());
   }


   m_OutTree->Fill();
   if(m_is_Zdc_Calib && ((ievent%10000)==0) ) m_OutTree->Write(); 
   ATH_MSG_DEBUG("Filling RunNumber="<<RunNumber);
   return StatusCode::SUCCESS;
}





StatusCode TrigRates::finalize(){
  return StatusCode::SUCCESS;
}





StatusCode TrigRates::CleaningCuts(bool &pass_cleaning_cuts){
   ATH_MSG_DEBUG("Starting CleaningCuts()");
   pass_cleaning_cuts=true;
   static int event_num=-1;
   event_num++;

//Event Info
//------------------------------------------------------  
//Must Pass good lumi-block
   const xAOD::EventInfo* l_EventInfo = nullptr;
   if(evtStore()->retrieve(l_EventInfo,m_EventInfo_key).isFailure()){
      ATH_MSG_ERROR(" Could not retrieve EventInfo with key "<<m_EventInfo_key.c_str());
      return StatusCode::FAILURE;
   }

   bool isMC=false;
   if(l_EventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
     isMC = true; // can do something with this later
     ATH_MSG_DEBUG("Event is MC");
   }  

   if(m_use_GRL && !isMC){
     if(!m_grlTool->passRunLB(*l_EventInfo)) {
       pass_cleaning_cuts=false;
       ATH_MSG_DEBUG("Event FAILED GRL");
       //std::cout<<"Event Failed GRL for run"<<l_EventInfo->runNumber()<<std::endl;
       return StatusCode::SUCCESS;
     }
     else{
       ATH_MSG_DEBUG("Event PASSED GRL");
       //std::cout<<"Event PASSED GRL for run"<<l_EventInfo->runNumber()<<std::endl;
     }
   }
   else {ATH_MSG_DEBUG("Not using GRL");}
//------------------------------------------------------  



//------------------------------------------------------  
//Event by Event Cleaning
    if(!isMC){
      int flag=0;
      if(l_EventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error  ) flag =1;
      if(l_EventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) flag+=2;
      if(l_EventInfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error  ) flag+=4;
      if(l_EventInfo->isEventFlagBitSet(xAOD::EventInfo::Core,18)               ) flag+=8;
      if(flag!=0){
          pass_cleaning_cuts=false;
          ATH_MSG_DEBUG("Event FAILED EventLevelCleaning :flag="<<flag);
          return StatusCode::SUCCESS;
      }
      else ATH_MSG_DEBUG("Event PASSED EventLevelCleaning");
    }
//------------------------------------------------------  



//Vertex
//------------------------------------------------------  
//Must have reconstructed vertex with |Z|<m_MaxZvtx
   if( (m_MaxZvtx>0) && (m_is_Zdc_Calib==false) ){
     const xAOD::VertexContainer *l_VertexContainer = nullptr;
     if(evtStore()->retrieve(l_VertexContainer,m_vtx_container_key).isFailure()){
        ATH_MSG_ERROR("Could not retrieve VxContainer with key "<<m_vtx_container_key.c_str());
        return StatusCode::FAILURE;
     }
     ATH_MSG_DEBUG("Num_Vertices="<<l_VertexContainer->size());
     if(l_VertexContainer->size()<=1){
       ATH_MSG_DEBUG("Event FAILED numvtx>=2 requirement");
       pass_cleaning_cuts=false;
     }

     int count=0;
     for(const auto* vtx : *l_VertexContainer){
        float z_vtx   =vtx->z();

        ATH_MSG_DEBUG("     z["<<count<<"]="<<z_vtx);
        if(count==0 && fabs(z_vtx)>m_MaxZvtx) {pass_cleaning_cuts=false;}
        count++;
     }
   }
//------------------------------------------------------  
   
   return StatusCode::SUCCESS;
}




//Triggers
void TrigRates::AddTriggerBranches(){
   m_ListOfTriggers      =m_trigTool->getChainGroup(m_Trigger_Chains       )->getListOfTriggers();
   m_ListOfMuonTriggers  =m_trigTool->getChainGroup(m_Muon_Trigger_Chains  )->getListOfTriggers();
   m_ListOfDiMuonTriggers=m_trigTool->getChainGroup(m_DiMuon_Trigger_Chains)->getListOfTriggers();

   ATH_MSG_INFO("STARTING LISTING TRIGGERS ::"<<m_Trigger_Chains);
   for(auto &trig : m_ListOfTriggers){
     ATH_MSG_INFO("Found Trigger : "<<trig.c_str());
   }
   for(auto &trig : m_ListOfMuonTriggers){
     ATH_MSG_INFO("Found Muon Trigger : "<<trig.c_str());
   }
   for(auto &trig : m_ListOfDiMuonTriggers){
     ATH_MSG_INFO("Found DiMuon Trigger : "<<trig.c_str());
   }
   ATH_MSG_INFO("FINISHED LISTING TRIGGERS ::"<<m_Trigger_Chains);



   int itrig=0;
   auto TriggerList=m_ListOfTriggers;
   if(m_store_single_muon   ) TriggerList.insert(TriggerList.end(),m_ListOfMuonTriggers.begin()  ,m_ListOfMuonTriggers.end()  );
   if(m_store_acoplanar_muon) TriggerList.insert(TriggerList.end(),m_ListOfDiMuonTriggers.begin(),m_ListOfDiMuonTriggers.end());

   for(auto trigger_name: TriggerList){
      //Branch storing trigger decision
      std::string BranchName="b_";BranchName+=trigger_name;
      std::string BranchType=BranchName+"/O";
      m_OutTree->Branch(BranchName.c_str(), &m_trigger_Flag[itrig],BranchType.c_str());

      bool save =false;
      bool save1=false;
      if(m_store_L1) save1=true;
      //if(trigger_name.find("mu4")!=std::string::npos) save=true;

      //Branch storing trigger prescale
      std::string BranchName1="f_"+trigger_name+"_prescale";
      std::string BranchType1=BranchName1+"/F";
      if(save) m_OutTree->Branch(BranchName1.c_str(), &m_trigger_prescale_Value[itrig],BranchType1.c_str());

      //Branch storing if trigger was prescaled in this event
      BranchName1="b_"+trigger_name+"_isPrescaled";
      BranchType1=BranchName1+"/O";
      if(save) m_OutTree->Branch(BranchName1.c_str(), &m_trigger_isPrescaled_Flag[itrig],BranchType1.c_str());

      //Branch storing L1-trigger before prescale in this event
      BranchName1="b_"+trigger_name+"_L1TBP";
      BranchType1=BranchName1+"/O";
      if(save || save1) m_OutTree->Branch(BranchName1.c_str(), &m_trigger_L1TBP_Flag[itrig],BranchType1.c_str());

      //Branch storing L1-trigger after prescale in this event
      BranchName1="b_"+trigger_name+"_L1TAP";
      BranchType1=BranchName1+"/O";
      if(save || save1) m_OutTree->Branch(BranchName1.c_str(), &m_trigger_L1TAP_Flag[itrig],BranchType1.c_str());

      //Branch storing L1-trigger after veto in this event
      BranchName1="b_"+trigger_name+"_L1TAV";
      BranchType1=BranchName1+"/O";
      if(save || save1) m_OutTree->Branch(BranchName1.c_str(), &m_trigger_L1TAV_Flag[itrig],BranchType1.c_str());

      //Branch storing if trigger was rerun
      BranchName1="b_"+trigger_name+"_isRerun";
      BranchType1=BranchName1+"/O";
      if(save) m_OutTree->Branch(BranchName1.c_str(), &m_trigger_isRerun_Flag[itrig],BranchType1.c_str());

      //Branch storing trigger decision with rerun
      BranchName1="b_"+trigger_name+"_rerun_decision";
      BranchType1=BranchName1+"/O";
      if(save) m_OutTree->Branch(BranchName1.c_str(), &m_trigger_rerun_Flag[itrig],BranchType1.c_str());

      m_trigger_map            [trigger_name]=&m_trigger_Flag            [itrig];
      m_trigger_prescale_map   [trigger_name]=&m_trigger_prescale_Value  [itrig];
      m_trigger_isPrescaled_map[trigger_name]=&m_trigger_isPrescaled_Flag[itrig];
      m_trigger_L1TBP_map      [trigger_name]=&m_trigger_L1TBP_Flag      [itrig];
      m_trigger_L1TAP_map      [trigger_name]=&m_trigger_L1TAP_Flag      [itrig];
      m_trigger_L1TAV_map      [trigger_name]=&m_trigger_L1TAV_Flag      [itrig];
      m_trigger_isRerun_map    [trigger_name]=&m_trigger_isRerun_Flag    [itrig];
      m_trigger_rerun_map      [trigger_name]=&m_trigger_rerun_Flag      [itrig];

      itrig++;
      if(itrig>=1000) {ATH_MSG_ERROR("Trigger count exceeded ");throw std::exception();}
   }
   itrig=0;
   if(m_store_single_muon){
     for(auto trigger_name: m_ListOfMuonTriggers){
        //Branch storing trigger matching for offline objects
        std::string BranchName="b_";BranchName+=trigger_name;
        std::string BranchName2="muon_"+BranchName;
        m_OutTree->Branch(BranchName2.c_str(), &m_trigger_match_Flag[itrig]);

        //Branch storing trigger matching for offline objects (without requirement that event passed trigger)
        BranchName2="muon_"+BranchName+"_V2";
        m_OutTree->Branch(BranchName2.c_str(), &m_trigger_match_Flag_V2[itrig]);

        //Branch storing trigger matching for offline objects (without requirement that event passed trigger)
        BranchName2="muon_"+BranchName+"_V3";
        m_OutTree->Branch(BranchName2.c_str(), &m_trigger_match_Flag_V3[itrig]);

        m_trigger_match_map      [trigger_name]=&m_trigger_match_Flag      [itrig];
        m_trigger_match_map_V2   [trigger_name]=&m_trigger_match_Flag_V2   [itrig];
        m_trigger_match_map_V3   [trigger_name]=&m_trigger_match_Flag_V3   [itrig];

        itrig++;
     }
   }
   if(m_store_acoplanar_muon){
     for(auto trigger_name: m_ListOfDiMuonTriggers){
        #if defined(__ATHENA_21p2__) || defined(__ATHENA_24p2__)
        //Branch storing trigger matching for offline objects
        std::string BranchName="b_";BranchName+=trigger_name;
        std::string BranchName2="dimuon_"+BranchName;
        m_OutTree->Branch(BranchName2.c_str(), &m_trigger_match_Flag[itrig]);
        #endif

        m_trigger_match_map      [trigger_name]=&m_trigger_match_Flag      [itrig];

        itrig++;
     }
   }
   if(itrig>=1000) {ATH_MSG_ERROR("Trigger count exceeded ");throw std::exception();}



   ATH_MSG_INFO("STARTING LISTING ALL TRIGGERS ::");
   auto ListOfAllTriggers=m_trigTool->getChainGroup(".*")->getListOfTriggers();
   for(auto &trig : ListOfAllTriggers){
     ATH_MSG_INFO("Found Trigger : "<<trig.c_str());
   }
   ATH_MSG_INFO("FINISHED LISTING ALL TRIGGERS ::");

}

StatusCode TrigRates::ProcessTriggers(bool &pass_trigger_cuts){
   for(const auto& trig_chain:m_trigger_map){
     std::string trigger_name=trig_chain.first ; //The name of the trigger 
     bool       *trigger_flag=trig_chain.second; //The bool that gets writted to the outputTree
     *trigger_flag                           =false;
     *m_trigger_prescale_map   [trigger_name]=1.0;
     *m_trigger_isPrescaled_map[trigger_name]=false;
     *m_trigger_L1TBP_map      [trigger_name]=false;
     *m_trigger_L1TAP_map      [trigger_name]=false;
     *m_trigger_L1TAV_map      [trigger_name]=false;
     *m_trigger_isRerun_map    [trigger_name]=false;
     *m_trigger_rerun_map      [trigger_name]=false;

     if  (m_trigTool->getListOfTriggers(trigger_name).empty()) {ATH_MSG_WARNING("Trigger "<<trigger_name.c_str()<<" Is Not Configured");}
     else                                                      {ATH_MSG_DEBUG  ("Trigger "<<trigger_name.c_str()<<" Is Configured"    );}

     //*m_trigger_prescale_map[trigger_name]=m_trigTool->getChainGroup(trigger_name)->getPrescale(TrigDefs::fullChain);
     *m_trigger_prescale_map[trigger_name]=m_trigTool->getChainGroup(trigger_name)->getPrescale(TrigDefs::Physics);
     //TODO which of the lines above should be used for CalibStream?

     if(m_trigTool->isPassed(trigger_name)){
        ATH_MSG_DEBUG("   Passed Trigger "<<trigger_name.c_str());
        *trigger_flag     =true;
        pass_trigger_cuts =true;
     }
     else{
        ATH_MSG_DEBUG("   Failed Trigger "<<trigger_name.c_str());
        *trigger_flag=false;
     }


     //Check if trigger was prescaled 
     const unsigned int triggerbits = m_trigTool->isPassedBits(trigger_name);
     bool efprescale = triggerbits & TrigDefs::EF_prescaled;
     // What about the L1?
     bool tbp = triggerbits&TrigDefs::L1_isPassedBeforePrescale;
     bool tap = triggerbits&TrigDefs::L1_isPassedAfterPrescale;
     bool tav = triggerbits&TrigDefs::L1_isPassedAfterVeto;
     // L1 is prescaled if tbp and not tap
     //bool l1prescale=tbp && !tap;
     bool l1prescale=tbp && !tav;
     bool chainPrescale= efprescale || l1prescale;
     *m_trigger_isPrescaled_map[trigger_name]=chainPrescale;
     *m_trigger_L1TBP_map      [trigger_name]=tbp;
     *m_trigger_L1TAP_map      [trigger_name]=tap;
     *m_trigger_L1TAV_map      [trigger_name]=tav;

     //rerun decision
     bool resurrected=triggerbits & TrigDefs::EF_resurrected;
     if(chainPrescale && resurrected){
       *m_trigger_isRerun_map[trigger_name]=true;
       if(m_trigTool->isPassed(trigger_name,TrigDefs::allowResurrectedDecision|TrigDefs::requireDecision)){
         *m_trigger_rerun_map[trigger_name]=true;
       }
       //if(m_trigTool->isPassed(trigger_name,TrigDefs::allowResurrectedDecision|TrigDefs::eventAccepted)){}
     }
   }
   return StatusCode::SUCCESS;
}






//EventInfo
void TrigRates::InitEventInfo(TTree *l_OutTree){
   l_OutTree->Branch("RunNumber"  , &RunNumber        ,"RunNumber/i");
   l_OutTree->Branch("lbn"        , &lumi_block       ,"lbn/i");
   l_OutTree->Branch("bcid"       , &bunch_crossing_id,"bcid/i");
   l_OutTree->Branch("eventNumber", &eventNumber      ,"eventNumber/l");
   l_OutTree->Branch("ActIntPerXing", &ActIntPerXing  ,"ActIntPerXing/F");
   l_OutTree->Branch("AvgIntPerXing", &AvgIntPerXing  ,"AvgIntPerXing/F");
   if(m_store_EventInfo&2){
     l_OutTree->Branch("ZdcEtA"              ,&ZdcEtA                ,"ZdcEtA/F");
     l_OutTree->Branch("ZdcEtC"              ,&ZdcEtC                ,"ZdcEtC/F");
     l_OutTree->Branch("is_pileup"           ,&m_is_pileup           ,"is_pileup/O");
     l_OutTree->Branch("is_oo_pileup"        ,&m_is_oo_pileup        ,"is_oo_pileup/O");
     l_OutTree->Branch("numtrk_HITight_pt500",&m_numtrk_HITight_pt500,"m_numtrk_HITight_pt500/I");
   }
   if(m_store_EventInfo&4){
     l_OutTree->Branch("NumTrackNoCuts"   ,&m_NumTrackNoCuts   ,"NumTrackNoCuts/I"   );
     l_OutTree->Branch("NumTrackPPMinBias",&m_NumTrackPPMinBias,"NumTrackPPMinBias/I");
     l_OutTree->Branch("NumTrackHILoose"  ,&m_NumTrackHILoose  ,"NumTrackHILoose/I"  );
     l_OutTree->Branch("NumTrackHITight"  ,&m_NumTrackHITight  ,"NumTrackHITight/I"  );

     l_OutTree->Branch("FCalET_aux"    ,&m_FCalET    ,"FCalET_aux/F"   );
     l_OutTree->Branch("FCalETP_aux"   ,&m_FCalETP   ,"FCalETP_aux/F"   );
     l_OutTree->Branch("FCalETN_aux"   ,&m_FCalETN   ,"FCalETN_aux/F"   );
   }
}

StatusCode TrigRates::ProcessEventInfo(){
   const xAOD::EventInfo* l_EventInfo = nullptr;
   if(evtStore()->retrieve(l_EventInfo,m_EventInfo_key).isFailure()){
      ATH_MSG_ERROR(" Could not retrieve EventInfo with key "<<m_EventInfo_key.c_str());
      return StatusCode::FAILURE;
   }

   //Set Variables written to output Tree
   RunNumber        =l_EventInfo->runNumber();
   lumi_block       =l_EventInfo->lumiBlock();
   bunch_crossing_id=l_EventInfo->bcid();
   eventNumber      =l_EventInfo->eventNumber();
   AvgIntPerXing    =l_EventInfo->averageInteractionsPerCrossing();
   ActIntPerXing    =l_EventInfo->actualInteractionsPerCrossing ();


   //Num tracks passig cuts
   //Only works in my personal derivation
   //Also used below 
   if(m_store_EventInfo&4){
      m_NumTrackNoCuts   =l_EventInfo->auxdata<int>("NumTrackNoCuts"   );
      m_NumTrackPPMinBias=l_EventInfo->auxdata<int>("NumTrackPPMinBias");
      m_NumTrackHILoose  =l_EventInfo->auxdata<int>("NumTrackHILoose"  );
      m_NumTrackHITight  =l_EventInfo->auxdata<int>("NumTrackHITight"  );

      m_FCalET   =l_EventInfo->auxdata<float>("FCalET"   );
      m_FCalETP  =l_EventInfo->auxdata<float>("FCalETP"  );
      m_FCalETN  =l_EventInfo->auxdata<float>("FCalETN"  );
   }


   //ZDC and pileup Info
   //Wont work in my personal derivation as the track container is thinned
   if(m_store_EventInfo&2){
     //-------------------------------------------------------------------------------
     ZdcEtA = 0;
     ZdcEtC = 0;
     #ifndef __ATHENA_24p2__
     m_ZDCAnalysisTool->reprocessZdc(); // Run the re-processing of the ZDC to get this information.
     const xAOD::ZdcModuleContainer *zdcSums = NULL;
     if(evtStore()->retrieve(zdcSums,  "ZdcSums"+m_ZdcAuxSuffix).isFailure()){
       ATH_MSG_ERROR(" Could not retrieve ZdcSums with key "<<"ZdcSums"+m_ZdcAuxSuffix);
       return StatusCode::FAILURE;
     }
     for(const auto zdcSum : *zdcSums){
       if(zdcSum->side() > 0) // Positive side of the ZDC 
         ZdcEtC+= zdcSum->auxdataConst<float>("CalibEnergy")/1000.0;//TODO should This be SideC?
       else // Negative side of the ZDC
         ZdcEtA+= zdcSum->auxdataConst<float>("CalibEnergy")/1000.0;//TODO should This be SideA?
     }
     #endif
     //-------------------------------------------------------------------------------


     //-------------------------------------------------------------------------------
     const xAOD::HIEventShapeContainer *hiue=nullptr;
     CHECK(evtStore()->retrieve(hiue, "CaloSums"));
     if(m_store_EventInfo&4){
       #ifndef __ATHENA_24p2__
       m_is_pileup    = m_HIPileupTool->is_pileup   ( *hiue, *zdcSums         );
       m_is_oo_pileup = m_HIPileupTool->is_Outpileup( *hiue, m_NumTrackHITight);
       #endif
     }
     else{
       const xAOD::TrackParticleContainer* tracks=nullptr;
       CHECK(evtStore()->retrieve( tracks, "InDetTrackParticles" ));
       m_numtrk_HITight_pt500=0;
       for (auto Track: *tracks) {
         if(m_trkSelTool_HITight->accept(Track) && Track->pt()>500) m_numtrk_HITight_pt500++;
       } 
       #ifndef __ATHENA_24p2__
       m_is_pileup    = m_HIPileupTool->is_pileup   ( *hiue, *zdcSums  );
       m_is_oo_pileup = m_HIPileupTool->is_Outpileup( *hiue, m_numtrk_HITight_pt500);
       #endif
     }
     //-------------------------------------------------------------------------------

   }




   return StatusCode::SUCCESS;
}



  float t_ZdcAmp      [2];
  float t_ZdcAmpErr   [2];
  float t_ZdcEnergy   [2];
  float t_ZdcEnergyErr[2];
  float t_ZdcTime     [2];
  short t_ZdcStatus   [2];
  unsigned int t_ZdcModuleMask;
  float t_ZdcTrigEff  [2];
  unsigned short t_ZdcLucrodTriggerSideAmp[2];

  //float t_RpdSubAmp[2][4][4];
  float t_RpdSubAmpSum         [2];
  float t_xDetCentroid         [2];
  float t_yDetCentroid         [2];
  float t_xCentroid            [2];
  float t_yCentroid            [2];
  float t_xDetCentroidUnsub    [2];
  float t_yDetCentroidUnsub    [2];
  float t_xDetRowCentroid      [2][4]={{0}};
  float t_yDetColCentroid      [2][4]={{0}};
  float t_xDetRowCentroidStdev [2];
  float t_yDetColCentroidStdev [2];
  float t_reactionPlaneAngle   [2];
  float t_cosDeltaReactionPlaneAngle;
  unsigned int t_centroidStatus[2];


void TrigRates::InitZdc(TTree *l_OutTree){
  #if defined(__ATHENA_24p2__)
  l_OutTree->Branch("zdc_ZdcAmp"       , &t_ZdcAmp      , "zdc_ZdcAmp[2]/F");
  l_OutTree->Branch("zdc_ZdcAmpErr"    , &t_ZdcAmpErr   , "zdc_ZdcAmpErr[2]/F");
  l_OutTree->Branch("zdc_ZdcEnergy"    , &t_ZdcEnergy   , "zdc_ZdcEnergy[2]/F");
  l_OutTree->Branch("zdc_ZdcEnergyErr" , &t_ZdcEnergyErr, "zdc_ZdcEnergyErr[2]/F");
  l_OutTree->Branch("zdc_ZdcTime"      , &t_ZdcTime     , "zdc_ZdcTime[2]/F");
  l_OutTree->Branch("zdc_ZdcStatus"    , &t_ZdcStatus   , "zdc_ZdcStatus[2]/S");
  //l_OutTree->Branch("zdc_ZdcTrigEff"   , &t_ZdcTrigEff  , "zdc_ZdcTrigEff[2]/F");
  l_OutTree->Branch("zdc_ZdcModuleMask", &t_ZdcModuleMask, "zdc_ZdcModuleMask/i");
  //l_OutTree->Branch("zdc_ZdcLucrodTriggerSideAmp",&t_ZdcLucrodTriggerSideAmp,"zdc_ZdcLucrodTriggerSideAmp[2]/S");

  if(m_store_Zdc & 2){
    //l_OutTree->Branch("zdc_RpdSubAmp"                 ,&t_RpdSubAmp                 ,"zdc_RpdSubAmp[2][4][4]/F");
    l_OutTree->Branch("zdc_RpdSubAmpSum"              ,&t_RpdSubAmpSum              ,"zdc_RpdSubAmpSum[2]/F");
    l_OutTree->Branch("zdc_xDetCentroid"              ,&t_xDetCentroid              ,"zdc_xDetCentroid[2]/F");
    l_OutTree->Branch("zdc_yDetCentroid"              ,&t_yDetCentroid              ,"zdc_yDetCentroid[2]/F");
    l_OutTree->Branch("zdc_xCentroid"                 ,&t_xCentroid                 ,"zdc_xCentroid[2]/F");
    l_OutTree->Branch("zdc_yCentroid"                 ,&t_yCentroid                 ,"zdc_yCentroid[2]/F");
    l_OutTree->Branch("zdc_xDetCentroidUnsub"         ,&t_xDetCentroidUnsub         ,"zdc_xDetCentroidUnsub[2]/F");
    l_OutTree->Branch("zdc_yDetCentroidUnsub"         ,&t_yDetCentroidUnsub         ,"zdc_yDetCentroidUnsub[2]/F");
    //l_OutTree->Branch("zdc_xDetRowCentroid"           ,&t_xDetRowCentroid           ,"zdc_xDetRowCentroid[2][4]/F");
    //l_OutTree->Branch("zdc_yDetColCentroid"           ,&t_yDetColCentroid           ,"zdc_yDetColCentroid[2][4]/F");
    l_OutTree->Branch("zdc_xDetRowCentroidStdev"      ,&t_xDetRowCentroidStdev      ,"zdc_xDetRowCentroidStdev[2]/F");
    l_OutTree->Branch("zdc_yDetColCentroidStdev"      ,&t_yDetColCentroidStdev      ,"zdc_yDetColCentroidStdev[2]/F");
    l_OutTree->Branch("zdc_reactionPlaneAngle"        ,&t_reactionPlaneAngle        ,"zdc_reactionPlaneAngle[2]/F");
    l_OutTree->Branch("zdc_cosDeltaReactionPlaneAngle",&t_cosDeltaReactionPlaneAngle,"zdc_cosDeltaReactionPlaneAngle/F");
    l_OutTree->Branch("zdc_centroidStatus"            ,&t_centroidStatus            ,"zdc_centroidStatus[2]/i");
  }
  #endif
}

StatusCode TrigRates::ProcessZdc(){
  #if defined(__ATHENA_24p2__)
  t_ZdcModuleMask =0;

  //m_ZDCAnalysisTool->reprocessZdc(); // Run the re-processing of the ZDC to get this information.

  const xAOD::ZdcModuleContainer *zdcSums = NULL;
  if(evtStore()->retrieve(zdcSums,  "ZdcSums").isFailure()){
    ATH_MSG_ERROR(" Could not retrieve ZdcSums with key "<<"ZdcSums");
    return StatusCode::FAILURE;
  }

  std::string auxSuffix=m_ZdcAuxSuffix;
  for(const auto zdcSum : *zdcSums){
    /*
    std::cout<<" side="         <<zdcSum->zdcSide()
             <<" Energy="      <<zdcSum->auxdataConst<float>("CalibEnergy"   + auxSuffix)
             <<" RpdSubAmpSum="<<zdcSum->auxdataConst<float>("RpdSubAmpSum"  + auxSuffix);
    std::cout<<std::endl;
    */

    int zdcside=zdcSum->zdcSide();
    //side==0 only for this varaible
    if(zdcside==0){
      t_cosDeltaReactionPlaneAngle=zdcSum->auxdataConst<float>("cosDeltaReactionPlaneAngle" + auxSuffix);
      continue;
    }

    const int iside=(zdcSum->zdcSide()>0)? 1:0;

    t_ZdcEnergy    [iside] = zdcSum->auxdataConst <float       >("CalibEnergy"   +auxSuffix );
    t_ZdcEnergyErr [iside] = zdcSum->auxdataConst <float       >("CalibEnergyErr"+auxSuffix );
    t_ZdcAmp       [iside] = zdcSum->auxdataConst <float       >("UncalibSum"    +auxSuffix );
    t_ZdcAmpErr    [iside] = zdcSum->auxdataConst <float       >("UncalibSumErr" +auxSuffix );
    t_ZdcTime      [iside] = zdcSum->auxdataConst <float       >("AverageTime"   +auxSuffix );
    t_ZdcStatus    [iside] = zdcSum->auxdataConst <unsigned int>("Status"        +auxSuffix );
    t_ZdcModuleMask       += (zdcSum->auxdataConst<unsigned int>("ModuleMask"    +auxSuffix ) << 4 * iside);

    if(m_store_Zdc & 2){
      t_RpdSubAmpSum     [iside] = zdcSum->auxdataConst<float>("RpdSubAmpSum"      + auxSuffix);
      t_xDetCentroid     [iside] = zdcSum->auxdataConst<float>("xDetCentroid"      + auxSuffix);
      t_yDetCentroid     [iside] = zdcSum->auxdataConst<float>("yDetCentroid"      + auxSuffix);
      t_xCentroid        [iside] = zdcSum->auxdataConst<float>("xCentroid"         + auxSuffix);
      t_yCentroid        [iside] = zdcSum->auxdataConst<float>("yCentroid"         + auxSuffix);
      t_xDetCentroidUnsub[iside] = zdcSum->auxdataConst<float>("xDetCentroidUnsub" + auxSuffix);
      t_yDetCentroidUnsub[iside] = zdcSum->auxdataConst<float>("yDetCentroidUnsub" + auxSuffix);

      std::vector<float> temp1=zdcSum->auxdataConst<std::vector<float>>("xDetRowCentroid" + auxSuffix);
      std::vector<float> temp2=zdcSum->auxdataConst<std::vector<float>>("yDetColCentroid" + auxSuffix);
      for (int row = 0; row < 4; row++) {
        t_xDetRowCentroid[iside][row] = temp1.at(row);
      }
      for (int col = 0; col < 4; col++) {
        t_yDetColCentroid[iside][col] = temp2.at(col);
      }
      t_xDetRowCentroidStdev[iside] = zdcSum->auxdataConst<float       >("xDetRowCentroidStdev" + auxSuffix);
      t_yDetColCentroidStdev[iside] = zdcSum->auxdataConst<float       >("yDetColCentroidStdev" + auxSuffix);
      t_reactionPlaneAngle  [iside] = zdcSum->auxdataConst<float       >("reactionPlaneAngle"   + auxSuffix);
      t_centroidStatus      [iside] = zdcSum->auxdataConst<unsigned int>("centroidStatus"       + auxSuffix);

    }
  }


  return StatusCode::SUCCESS;
  #else
  ATH_MSG_ERROR("Zdc analysis is not implemented for this release");
  return StatusCode::FAILURE;
  #endif
}


//Vertex
void TrigRates::InitVertex(TTree *l_OutTree){
     l_OutTree->Branch("vtx_z"     ,&m_vtx_z     );
     l_OutTree->Branch("vtx_x"     ,&m_vtx_x     );
     l_OutTree->Branch("vtx_y"     ,&m_vtx_y     );
     l_OutTree->Branch("vtx_ntrk"  ,&m_vtx_ntrk  );
}

StatusCode TrigRates::ProcessVertex(){
   const xAOD::VertexContainer *l_VertexContainer = nullptr;
   if(evtStore()->retrieve(l_VertexContainer,m_vtx_container_key).isFailure()){
      ATH_MSG_ERROR("Could not retrieve VxContainer with key "<<m_vtx_container_key.c_str());
      return StatusCode::FAILURE;
   }

   m_vtx_z     .clear();
   m_vtx_x     .clear();
   m_vtx_y     .clear();
   m_vtx_ntrk  .clear();
   for(const auto* vtx : *l_VertexContainer){
      m_vtx_z     .push_back(vtx->z());
      m_vtx_x     .push_back(vtx->x());
      m_vtx_y     .push_back(vtx->y());
      m_vtx_ntrk  .push_back(vtx->nTrackParticles());
    }
   return StatusCode::SUCCESS;
}


//Muons
void TrigRates::InitMuons(TTree *l_OutTree){
  if(m_store_single_muon){
    #if defined(__ATHENA_24p2__)
    l_OutTree->Branch("muon_match_mu4roi",&m_muon_match_mu4roi);
    l_OutTree->Branch("muon_match_mu6roi",&m_muon_match_mu6roi);
    #endif

    l_OutTree->Branch("muon_pt"          ,&m_muon_pt          );
    l_OutTree->Branch("muon_eta"         ,&m_muon_eta         );
    l_OutTree->Branch("muon_phi"         ,&m_muon_phi         );

    l_OutTree->Branch("muon_pt_precorr"  ,&m_muon_pt_precorr  );
    l_OutTree->Branch("muon_eta_precorr" ,&m_muon_eta_precorr );
    l_OutTree->Branch("muon_phi_precorr" ,&m_muon_phi_precorr );

    l_OutTree->Branch("muon_quality"     ,&m_muon_qual        );
    l_OutTree->Branch("muon_deltaP_overP",&m_muon_deltaP_overP);

    l_OutTree->Branch("muon_eff_corr_medium_data",&m_muon_eff_corr_medium_data);
    l_OutTree->Branch("muon_eff_corr_medium_mc"  ,&m_muon_eff_corr_medium_mc  );
    l_OutTree->Branch("muon_eff_SF_medium"       ,&m_muon_eff_SF_medium       );
    l_OutTree->Branch("muon_eff_corr_tight_data" ,&m_muon_eff_corr_tight_data );
    l_OutTree->Branch("muon_eff_corr_tight_mc"   ,&m_muon_eff_corr_tight_mc   );
    l_OutTree->Branch("muon_eff_SF_tight"        ,&m_muon_eff_SF_tight        );

    l_OutTree->Branch("muon_d0"          ,&m_muon_d0          );
    l_OutTree->Branch("muon_d0_err"      ,&m_muon_d0_err      );
    l_OutTree->Branch("muon_z0"          ,&m_muon_z0          );
    l_OutTree->Branch("muon_trk_p"       ,&m_muon_trk_p       );
    l_OutTree->Branch("muon_trk_pt"      ,&m_muon_trk_pt      );
    l_OutTree->Branch("muon_trk_eta"     ,&m_muon_trk_eta     );
    l_OutTree->Branch("muon_trk_phi"     ,&m_muon_trk_phi     );
    if(m_store_tracks>=Track::StoreBasic) l_OutTree->Branch("muon_trk_index"   ,&m_muon_trk_index   );

    l_OutTree->Branch("muon_me_p"        ,&m_muon_me_p        );

    l_OutTree->Branch("muon_elosstype"   ,&m_muon_elosstype   );
    l_OutTree->Branch("muon_mspt"        ,&m_muon_mspt        );
    l_OutTree->Branch("muon_msp"         ,&m_muon_msp         );
    l_OutTree->Branch("muon_ms_phi_hits" ,&m_muon_ms_phi_hits );
    l_OutTree->Branch("muon_ms_eta_hits" ,&m_muon_ms_eta_hits );

    if(m_use_trigger)  l_OutTree->Branch("muon_trig_match"  ,&m_muon_trig_match  );

    if(m_store_truth){
      l_OutTree->Branch("muon_truth_index"    ,&m_muon_truth_index    ); 
    }
    if(m_store_muon_truth){
      l_OutTree->Branch("muon_truth_prob"     ,&m_muon_truth_prob     ); 
      l_OutTree->Branch("muon_truth_barcode"  ,&m_muon_truth_barcode  ); 
      l_OutTree->Branch("muon_truth_IsPrimary",&m_muon_truth_IsPrimary); 
      l_OutTree->Branch("muon_truth_pt"       ,&m_muon_truth_pt       );
      l_OutTree->Branch("muon_truth_eta"      ,&m_muon_truth_eta      );
      l_OutTree->Branch("muon_truth_phi"      ,&m_muon_truth_phi      );
      l_OutTree->Branch("muon_truth_charge"   ,&m_muon_truth_charge   );
      l_OutTree->Branch("muon_truth_id"       ,&m_muon_truth_id       );
      l_OutTree->Branch("muon_truth_quality"  ,&m_muon_truth_quality  );
      l_OutTree->Branch("muon_truth_status"   ,&m_muon_truth_status   );
    }
  }

  if(m_store_acoplanar_muon){
    l_OutTree->Branch("muon_pair_acop"              ,&m_muon_pair_acop              );
    l_OutTree->Branch("muon_pair_d0"                ,&m_muon_pair_d0                );
   
    l_OutTree->Branch("muon_pair_muon1_index"       ,&m_muon_pair_muon1_index       );
    l_OutTree->Branch("muon_pair_muon2_index"       ,&m_muon_pair_muon2_index       );

    l_OutTree->Branch("muon_pair_muon1_pt"          ,&m_muon_pair_muon1_pt          );
    l_OutTree->Branch("muon_pair_muon1_eta"         ,&m_muon_pair_muon1_eta         );
    l_OutTree->Branch("muon_pair_muon1_phi"         ,&m_muon_pair_muon1_phi         );
    //l_OutTree->Branch("muon_pair_muon1_pt_precorr"  ,&m_muon_pair_muon1_pt_precorr  );
    //l_OutTree->Branch("muon_pair_muon1_eta_precorr" ,&m_muon_pair_muon1_eta_precorr );
    //l_OutTree->Branch("muon_pair_muon1_phi_precorr" ,&m_muon_pair_muon1_phi_precorr );
    l_OutTree->Branch("muon_pair_muon1_quality"     ,&m_muon_pair_muon1_qual        );
    //l_OutTree->Branch("muon_pair_muon1_deltaP_overP",&m_muon_pair_muon1_deltaP_overP);
    l_OutTree->Branch("muon_pair_muon1_d0"          ,&m_muon_pair_muon1_d0          );
    l_OutTree->Branch("muon_pair_muon1_d0_err"      ,&m_muon_pair_muon1_d0_err      );
    l_OutTree->Branch("muon_pair_muon1_z0"          ,&m_muon_pair_muon1_z0          );
    if(m_store_tracks>=Track::StoreBasic) l_OutTree->Branch("muon_pair_muon1_trk_index"   ,&m_muon_pair_muon1_trk_index   );
    l_OutTree->Branch("muon_pair_muon1_mspt"        ,&m_muon_pair_muon1_mspt        );
    l_OutTree->Branch("muon_pair_muon1_msp"         ,&m_muon_pair_muon1_msp         );
    //l_OutTree->Branch("muon_pair_muon1_elosstype"   ,&m_muon_pair_muon1_elosstype   );
    l_OutTree->Branch("muon_pair_muon1_trk_pt"      ,&m_muon_pair_muon1_trk_pt      );
    l_OutTree->Branch("muon_pair_muon1_trk_eta"     ,&m_muon_pair_muon1_trk_eta     );
    l_OutTree->Branch("muon_pair_muon1_trk_phi"     ,&m_muon_pair_muon1_trk_phi     );

    l_OutTree->Branch("muon_pair_muon2_pt"          ,&m_muon_pair_muon2_pt          );
    l_OutTree->Branch("muon_pair_muon2_eta"         ,&m_muon_pair_muon2_eta         );
    l_OutTree->Branch("muon_pair_muon2_phi"         ,&m_muon_pair_muon2_phi         );
    //l_OutTree->Branch("muon_pair_muon2_pt_precorr"  ,&m_muon_pair_muon2_pt_precorr  );
    //l_OutTree->Branch("muon_pair_muon2_eta_precorr" ,&m_muon_pair_muon2_eta_precorr );
    //l_OutTree->Branch("muon_pair_muon2_phi_precorr" ,&m_muon_pair_muon2_phi_precorr );
    l_OutTree->Branch("muon_pair_muon2_quality"     ,&m_muon_pair_muon2_qual        );
    //l_OutTree->Branch("muon_pair_muon2_deltaP_overP",&m_muon_pair_muon2_deltaP_overP);
    l_OutTree->Branch("muon_pair_muon2_d0"          ,&m_muon_pair_muon2_d0          );
    l_OutTree->Branch("muon_pair_muon2_d0_err"      ,&m_muon_pair_muon2_d0_err      );
    l_OutTree->Branch("muon_pair_muon2_z0"          ,&m_muon_pair_muon2_z0          );
    if(m_store_tracks>=Track::StoreBasic) l_OutTree->Branch("muon_pair_muon2_trk_index"   ,&m_muon_pair_muon2_trk_index   );
    l_OutTree->Branch("muon_pair_muon2_mspt"        ,&m_muon_pair_muon2_mspt        );
    l_OutTree->Branch("muon_pair_muon2_msp"         ,&m_muon_pair_muon2_msp         );
    //l_OutTree->Branch("muon_pair_muon2_elosstype"   ,&m_muon_pair_muon2_elosstype   );
    l_OutTree->Branch("muon_pair_muon2_trk_pt"      ,&m_muon_pair_muon2_trk_pt      );
    l_OutTree->Branch("muon_pair_muon2_trk_eta"     ,&m_muon_pair_muon2_trk_eta     );
    l_OutTree->Branch("muon_pair_muon2_trk_phi"     ,&m_muon_pair_muon2_trk_phi     );
  }
}
void TrigRates::ClearMuons()
{
   #if defined(__ATHENA_24p2__)
   m_muon_match_mu4roi.clear();
   m_muon_match_mu6roi.clear();
   #endif

   m_muon_pt          .clear();
   m_muon_eta         .clear();
   m_muon_phi         .clear();

   m_muon_pt_precorr  .clear();
   m_muon_eta_precorr .clear();
   m_muon_phi_precorr .clear();

   m_muon_qual        .clear();
   m_muon_deltaP_overP.clear();

   m_muon_eff_corr_medium_data.clear();
   m_muon_eff_corr_medium_mc  .clear();
   m_muon_eff_SF_medium       .clear();
   m_muon_eff_corr_tight_data .clear();
   m_muon_eff_corr_tight_mc   .clear();
   m_muon_eff_SF_tight        .clear();

   m_muon_d0          .clear();
   m_muon_d0_err      .clear();
   m_muon_z0          .clear();
   m_muon_trk_p       .clear();
   m_muon_trk_pt      .clear();
   m_muon_trk_eta     .clear();
   m_muon_trk_phi     .clear();
   m_muon_trk_index   .clear();

   m_muon_me_p        .clear();

   m_muon_elosstype   .clear();
   m_muon_mspt        .clear();
   m_muon_msp         .clear();
   m_muon_ms_phi_hits .clear();
   m_muon_ms_eta_hits .clear();

   m_muon_trig_match  .clear();
   if(m_store_truth){
     m_muon_truth_index_temp1.clear();
     m_muon_truth_index      .clear();
   }
   if(m_store_muon_truth){
     m_muon_truth_prob       .clear();
     m_muon_truth_barcode    .clear();
     m_muon_truth_IsPrimary  .clear();
     m_muon_truth_pt         .clear();
     m_muon_truth_eta        .clear();
     m_muon_truth_phi        .clear();
     m_muon_truth_charge     .clear();
     m_muon_truth_id         .clear();
     m_muon_truth_quality    .clear();
     m_muon_truth_status     .clear();
   }
   if(m_store_acoplanar_muon){
     m_muon_pair_acop              .clear();
     m_muon_pair_d0                .clear();

     m_muon_pair_muon1_index       .clear();
     m_muon_pair_muon2_index       .clear();

     m_muon_pair_muon1_pt          .clear();
     m_muon_pair_muon1_eta         .clear();
     m_muon_pair_muon1_phi         .clear();
     m_muon_pair_muon1_pt_precorr  .clear();
     m_muon_pair_muon1_eta_precorr .clear();
     m_muon_pair_muon1_phi_precorr .clear();
     m_muon_pair_muon1_qual        .clear();
     m_muon_pair_muon1_deltaP_overP.clear();
     m_muon_pair_muon1_d0          .clear();
     m_muon_pair_muon1_d0_err      .clear();
     m_muon_pair_muon1_z0          .clear();
     m_muon_pair_muon1_trk_index   .clear();
     m_muon_pair_muon1_mspt        .clear();
     m_muon_pair_muon1_msp         .clear();
     m_muon_pair_muon1_elosstype   .clear();
     m_muon_pair_muon1_trk_pt      .clear();
     m_muon_pair_muon1_trk_eta     .clear();
     m_muon_pair_muon1_trk_phi     .clear();
     m_muon_pair_muon2_pt          .clear();
     m_muon_pair_muon2_eta         .clear();
     m_muon_pair_muon2_phi         .clear();
     m_muon_pair_muon2_pt_precorr  .clear();
     m_muon_pair_muon2_eta_precorr .clear();
     m_muon_pair_muon2_phi_precorr .clear();
     m_muon_pair_muon2_qual        .clear();
     m_muon_pair_muon2_deltaP_overP.clear();
     m_muon_pair_muon2_d0          .clear();
     m_muon_pair_muon2_d0_err      .clear();
     m_muon_pair_muon2_z0          .clear();
     m_muon_pair_muon2_trk_index   .clear();
     m_muon_pair_muon2_mspt        .clear();
     m_muon_pair_muon2_msp         .clear();
     m_muon_pair_muon2_elosstype   .clear();
     m_muon_pair_muon2_trk_pt      .clear();
     m_muon_pair_muon2_trk_eta     .clear();
     m_muon_pair_muon2_trk_phi     .clear();
   }
}

StatusCode TrigRates::ProcessMuons(){
   ClearMuons();

   const xAOD::MuonContainer *l_muons;
   CHECK(evtStore()->retrieve(l_muons,m_muons_key));

   const xAOD::VertexContainer *l_VertexContainer = nullptr;
   CHECK(evtStore()->retrieve(l_VertexContainer,m_vtx_container_key));
   const auto* priVtx = *(l_VertexContainer->cbegin());
   const float z_vtx = priVtx->z();


   #if defined(__ATHENA_24p2__)
   const xAOD::MuonContainer *l_muons_trig_roi;
   CHECK(evtStore()->retrieve(l_muons_trig_roi,m_hlt_muons_key));
   #endif


   std::map<int,const xAOD::IParticle*> muon_map;//stores link between muons and original muon in the Container
   int imuon=0;
   for(auto Muon1:*l_muons){
     const xAOD::Muon* Muon = Muon1;
     xAOD::Muon* muon = nullptr;

     float pt_precorr =Muon->pt ();
     float eta_precorr=Muon->eta();
     float phi_precorr=Muon->phi();

     //-----------------------------------------------------------------
     //muonCalibrationAndSmearingTool is not fully implemented
     if(m_ApplyMuonCalibrations){
       if(!m_muonCalibrationAndSmearingTool->correctedCopy(*Muon,muon)) {
         ATH_MSG_INFO ("execute(): Problem with Muon Calibration And Smearing Tool (Error or OutOfValidityRange) ");
         if(muon){delete muon; muon=nullptr;}//in this case, continue with original muon from container
       }
       else{
         Muon=muon;
       }
     }
     //-----------------------------------------------------------------


     //-----------------------------------------------------------------
     int quality=0;
     xAOD::Muon::Quality my_quality = m_muonSelection->getQuality(*Muon);
     if(Muon->muonType()==xAOD::Muon::MuonType::Combined)  quality+= 1;//THERE WAS BUG HERE FIXED November 7 2018
     if(my_quality<=xAOD::Muon::VeryLoose)                 quality+= 2;
     if(my_quality<=xAOD::Muon::Loose)                     quality+= 4;
     if(my_quality<=xAOD::Muon::Medium)                    quality+= 8;
     if(my_quality<=xAOD::Muon::Tight)                     quality+=16;
     if(m_muonSelection->passedIDCuts(*Muon))              quality+=32;
     const xAOD::TrackParticle*idTrk= Muon->trackParticle(xAOD::Muon::InnerDetectorTrackParticle);//ID Track
     const xAOD::TrackParticle*meTrk= Muon->trackParticle(xAOD::Muon::ExtrapolatedMuonSpectrometerTrackParticle);//ME Track
     if(meTrk) quality+=64;
     if(idTrk) quality+=128;
     if(m_muonSelection->passedMuonCuts(*Muon))            quality+=256;
     //-----------------------------------------------------------------


     //-----------------------------------------------------------------
     float eff_corr_medium_data=0,eff_corr_medium_mc=0;
     float eff_corr_tight_data =0,eff_corr_tight_mc =0;
     float eff_SF_medium       =0,eff_SF_tight      =0;
     CP::CorrectionCode cr1=CP::CorrectionCode::ErrorCode::Ok,
                        cr2=CP::CorrectionCode::ErrorCode::Ok,
                        cr3=CP::CorrectionCode::ErrorCode::Ok;
     if(my_quality<=xAOD::Muon::Medium && m_use_effi_SF_tool){
       cr1=m_effi_SF_tool_medium->getEfficiencyScaleFactor(*Muon,eff_SF_medium       );
       cr2=m_effi_SF_tool_medium->getDataEfficiency       (*Muon,eff_corr_medium_data);
       cr3=m_effi_SF_tool_medium->getMCEfficiency         (*Muon,eff_corr_medium_mc  );
     }
     if(my_quality<=xAOD::Muon::Tight && m_use_effi_SF_tool){
       cr1=m_effi_SF_tool_tight ->getEfficiencyScaleFactor(*Muon,eff_SF_tight        );
       cr2=m_effi_SF_tool_tight ->getDataEfficiency       (*Muon,eff_corr_tight_data );
       cr3=m_effi_SF_tool_tight ->getMCEfficiency         (*Muon,eff_corr_tight_mc   );
     }
     //TODO Throw exception here?
     if(cr1!=CP::CorrectionCode::ErrorCode::Ok || 
        cr2!=CP::CorrectionCode::ErrorCode::Ok ||
        cr3!=CP::CorrectionCode::ErrorCode::Ok){
       eff_corr_medium_data=0,eff_corr_medium_mc=0;
       eff_corr_tight_data =0,eff_corr_tight_mc =0;
       eff_SF_medium       =0,eff_SF_tight      =0;
     }
     //-----------------------------------------------------------------


     //-----------------------------------------------------------------
     float deltaP_overP=1000;
     if(!meTrk || !idTrk){
       if(!meTrk && (quality&1)) {ATH_MSG_ERROR("meTrk Not Found for Combined muon");}
       if(!idTrk && (quality&1)) {ATH_MSG_ERROR("idTrk Not Found for Combined muon");}
     }
     else{
       deltaP_overP=(idTrk->p4().P() - meTrk->p4().P())/idTrk->p4().P();
     }
     //-----------------------------------------------------------------

     
     //-----------------------------------------------------------------
     float d0=1000,z0=1000,muon_trk_p=-1000,muon_me_p=-1000,muon_trk_pt=0,muon_trk_eta=1000,muon_trk_phi=0;
     int trk_index=-1;
     float d0_err=0;
     if(idTrk){
       xAOD::ParametersCovMatrix_t covmat=idTrk->definingParametersCovMatrix();
       d0          =idTrk->d0();
       d0_err      =sqrt(fabs(covmat(0,0)));
       z0          =idTrk->z0()+idTrk->vz() - z_vtx;

       muon_trk_p  =idTrk->p4().P();
       muon_trk_pt =idTrk->pt ();
       muon_trk_eta=idTrk->eta();
       muon_trk_phi=idTrk->phi();
       if(m_store_tracks>=Track::StoreBasic){
         if(m_track_index_temp.find(idTrk) != m_track_index_temp.end()) trk_index=m_track_index_temp[idTrk];
       }
     }
     if(meTrk)muon_me_p   =meTrk->p4().P();
     //-----------------------------------------------------------------


     //-----------------------------------------------------------------
     //see http://acode-browser2.usatlas.bnl.gov/lxr/source/r21/atlas/Event/xAOD/xAODMuon/
     //           xAODMuon/versions/Muon_v1.h (or latest file)
     ///// Defines how the energy loss was handled for this muon 
     //enum EnergyLossType { 
     //   Parametrized=0,  
     //   NotIsolated=1,  //!< Reconstruction configured to use the parametrization w/o looking in the calo (eg calo off)
     //   MOP=2,          //!< Measurement found to be compatible with most probable value --> 
     //                   //   mop used as more reliable at this region of the eloss
     //   Tail=3,         //!< Measured eloss significantly higher than mop --> the calo measurement used
     //   FSRcandidate=4  //!< In standalone reconstruction the Tail option was used. 
     //                   //but an imbalance is observed when comparing Pstandalone and Pinnerdetector (Pstandalone>Pinnerdetector) 
     //                   //--> if using the mop resolves the imbalance the excess energy loss is stored as 
     //                   //    fsrEnergy and the mop is used as the eloss.
     //};
     int elosstype=0;//Muon->auxdata<uint8_t>("energyLossType");
     //-----------------------------------------------------------------


     //-----------------------------------------------------------------
     float ms_pt=-1000,ms_p=-1000;
     uint8_t ms_phi_hits=0,ms_eta_hits=0;
     const xAOD::TrackParticle*msTrk= Muon->trackParticle(xAOD::Muon::MuonSpectrometerTrackParticle);//MS Track
     if(msTrk){
       ms_pt=msTrk->pt();
       ms_p =msTrk->p4().P();
       bool b1=msTrk->summaryValue(ms_phi_hits,xAOD::SummaryType::numberOfPhiLayers);
       bool b2=msTrk->summaryValue(ms_eta_hits,xAOD::SummaryType::numberOfTriggerEtaLayers);
       if(b1!=true || b2!=true){
         ATH_MSG_ERROR("Coudnot determine Phi/Eta Layer hits for muon");
       }
     }
     //-----------------------------------------------------------------


     //-----------------------------------------------------------------
     bool is_matched=false;
     if(m_store_single_muon && m_use_trigger){
       for(std::string trig_chain_name:m_ListOfMuonTriggers){
         #if defined(__ATHENA_21p2__) || defined(__ATHENA_24p2__)
         if(m_matchTool->match(*Muon,trig_chain_name,0.1)==true){ //}
         #else
         if(m_matchTool->match(Muon,trig_chain_name,0.1)==true){
         #endif
           m_trigger_match_map[trig_chain_name]->push_back(true);
           is_matched=true; 
         }
         else{
           //std::cout<<trig_chain_name<<std::endl;
           m_trigger_match_map[trig_chain_name]->push_back(false);
         }

         //https://gitlab.cern.ch/atlas/athena/-/blob/21.2/Trigger/TrigAnalysis/TriggerMatchingTool/Root/MatchingTool.cxx#L146
         m_trigger_match_map_V3[trig_chain_name]->push_back(m_matchTool->match(*Muon,trig_chain_name,0.1,true));
       }
     }
     //-----------------------------------------------------------------


     //-----------------------------------------------------------------
     //Manual matching,  this is done as the normal tool requires trigger to pass "Physics" and cannot be used for rerun triggers
     //https://gitlab.cern.ch/atlas/athena/blob/21.0/Trigger/TrigAnalysis/TrigMuonMatching/Root/TrigMuonMatching.cxx#L266
     if(m_store_single_muon && m_use_trigger){
       for(std::string trig_chain_name:m_ListOfMuonTriggers){
         bool match=false;
         #ifndef __ATHENA_24p2__ //This part fails for Rel>=22 currently (works for Rel21)
         auto cg = m_trigTool->getChainGroup(trig_chain_name);     
         auto fc = cg->features();
         auto MuFeatureContainers=fc.containerFeature<xAOD::MuonContainer>();
         //const std::vector< Trig::Feature<xAOD::MuonContainer> > MuFeatureContainers = fc.get<xAOD::MuonContainer>();
         for(auto mucont : MuFeatureContainers){
           for(auto mu : *mucont.cptr()){
             double deta=fabs(Muon->eta()-mu->eta());
             double dphi=fabs(Muon->phi()-mu->phi());
             if(dphi> TMath::Pi()){
               dphi=fabs(2*TMath::Pi()-dphi);
             }
             Double_t dr=sqrt(deta*deta + dphi*dphi);
             if(dr<0.1) {match=true;break;}
           }     
           if(match) break;
         }
         #endif
         m_trigger_match_map_V2[trig_chain_name]->push_back(match);
       }
     }
     //-----------------------------------------------------------------


     //-----------------------------------------------------------------
     int   truth_index=-1,truth_id=0,truth_quality=-1;
     float match_prob =0,truth_pt=0,truth_eta=1000,truth_phi=1000,truth_charge=1000;
     float barcode    =-1;
     bool  is_primary =false,truth_status=false;

     if(m_store_truth || m_store_muon_truth){
         if(idTrk){
           const                 ElementLink<xAOD::TruthParticleContainer> ptruthContainer=
                 (idTrk->auxdata<ElementLink<xAOD::TruthParticleContainer>>("truthParticleLink" ));
           if(ptruthContainer.isValid()){
             const xAOD::TruthParticle *associated_truth=*ptruthContainer;
             m_muon_truth_index_temp1   [associated_truth]=imuon; 

             truth_index=-1;
             match_prob =idTrk->auxdata<float>("truthMatchProbability"); 
             barcode    =associated_truth->barcode();
             is_primary =IsPrimaryParticle(associated_truth);

             //added 30/9/2021
             truth_pt    =associated_truth->pt    ();
             truth_eta   =(fabs(truth_pt)>0.00001)?associated_truth->eta():1000;
             truth_phi   =associated_truth->phi   ();
             truth_charge=associated_truth->charge();
             truth_id    =associated_truth->pdgId ();
             truth_quality =(associated_truth->isStrangeBaryon() ||associated_truth->barcode()<=0)?-1:1;
             truth_status  =(associated_truth->status()!=1)? false:true;
           }
         }
     }
     //-----------------------------------------------------------------


     //-----------------------------------------------------------------
     float pt=Muon->pt();
     if(Muon->charge()<0) pt=-pt;
     //-----------------------------------------------------------------

     #if defined(__ATHENA_24p2__)
      bool match_mu4roi=false;
      bool match_mu6roi=false;
      for(auto muon_trig:*l_muons_trig_roi){
        float pt_ =muon_trig->pt ();
        float eta_=muon_trig->eta();
        float phi_=muon_trig->phi();
        float dr2=pow(eta_precorr-eta_,2.0) + pow(phi_precorr-phi_,2.0);
        if(dr2<0.01){
          if(pt_>4000.0) match_mu4roi=true;
          if(pt_>6000.0) match_mu6roi=true;
        }
      }
      m_muon_match_mu4roi.push_back(match_mu4roi);
      m_muon_match_mu6roi.push_back(match_mu6roi);
     #endif


     m_muon_pt          .push_back(pt);
     m_muon_eta         .push_back(Muon->eta());
     m_muon_phi         .push_back(Muon->phi());

     m_muon_pt_precorr  .push_back(pt_precorr );
     m_muon_eta_precorr .push_back(eta_precorr);
     m_muon_phi_precorr .push_back(phi_precorr);

     m_muon_qual        .push_back(quality);
     m_muon_deltaP_overP.push_back(deltaP_overP);

     m_muon_eff_corr_medium_data.push_back(eff_corr_medium_data);
     m_muon_eff_corr_medium_mc  .push_back(eff_corr_medium_mc  );
     m_muon_eff_SF_medium       .push_back(eff_SF_medium       );
     m_muon_eff_corr_tight_data .push_back(eff_corr_tight_data );
     m_muon_eff_corr_tight_mc   .push_back(eff_corr_tight_mc   );
     m_muon_eff_SF_tight        .push_back(eff_SF_tight        );

     m_muon_d0          .push_back(d0);
     m_muon_d0_err      .push_back(d0_err);
     m_muon_z0          .push_back(z0);
     m_muon_trk_p       .push_back(muon_trk_p);
     m_muon_trk_pt      .push_back(muon_trk_pt);
     m_muon_trk_eta     .push_back(muon_trk_eta);
     m_muon_trk_phi     .push_back(muon_trk_phi);
     m_muon_trk_index   .push_back(trk_index);

     m_muon_me_p        .push_back(muon_me_p);

     m_muon_elosstype   .push_back(elosstype);
     m_muon_mspt        .push_back(ms_pt);
     m_muon_msp         .push_back(ms_p);
     m_muon_ms_phi_hits .push_back(ms_phi_hits);
     m_muon_ms_eta_hits .push_back(ms_eta_hits);

     m_muon_trig_match  .push_back(is_matched);

     if(m_store_truth){
       m_muon_truth_index    .push_back(truth_index  );
     }
     if(m_store_muon_truth){
       m_muon_truth_prob     .push_back(match_prob   ); 
       m_muon_truth_barcode  .push_back(barcode      );
       m_muon_truth_IsPrimary.push_back(is_primary   );
       m_muon_truth_pt       .push_back(truth_pt     );
       m_muon_truth_eta      .push_back(truth_eta    );
       m_muon_truth_phi      .push_back(truth_phi    );
       m_muon_truth_charge   .push_back(truth_charge );
       m_muon_truth_id       .push_back(truth_id     );
       m_muon_truth_quality  .push_back(truth_quality);
       m_muon_truth_status   .push_back(truth_status );
     }
     muon_map[imuon]=Muon1;

     if(muon){delete muon; muon = nullptr;}
     imuon++;
   }

   //Acoplanar or dimuon branches
   const double PI=acos(-1.0);
   if(m_store_acoplanar_muon){
     const int N=m_muon_pt.size();
     for(int i=0;i<N;i++){
       for(int j=i+1;j<N;j++){
         float phi1=m_muon_phi[i];
         float phi2=m_muon_phi[j];
         float dphi=phi1-phi2;
         float acop=1.0-fabs(atan2(sin(dphi),cos(dphi)))/PI;
         float d0_ =sqrt( pow(m_muon_d0[i],2.0) + pow(m_muon_d0[j],2.0) );

         m_muon_pair_acop              .push_back(acop);
         m_muon_pair_d0                .push_back(d0_);

         m_muon_pair_muon1_index       .push_back(i);
         m_muon_pair_muon2_index       .push_back(j);

         m_muon_pair_muon1_pt          .push_back(m_muon_pt          [i]);
         m_muon_pair_muon1_eta         .push_back(m_muon_eta         [i]);
         m_muon_pair_muon1_phi         .push_back(m_muon_phi         [i]);
         m_muon_pair_muon1_pt_precorr  .push_back(m_muon_pt_precorr  [i]);
         m_muon_pair_muon1_eta_precorr .push_back(m_muon_eta_precorr [i]);
         m_muon_pair_muon1_phi_precorr .push_back(m_muon_phi_precorr [i]);
         m_muon_pair_muon1_qual        .push_back(m_muon_qual        [i]);
         m_muon_pair_muon1_deltaP_overP.push_back(m_muon_deltaP_overP[i]);
         m_muon_pair_muon1_d0          .push_back(m_muon_d0          [i]);
         m_muon_pair_muon1_d0_err      .push_back(m_muon_d0_err      [i]);
         m_muon_pair_muon1_z0          .push_back(m_muon_z0          [i]);
         m_muon_pair_muon1_trk_index   .push_back(m_muon_trk_index   [i]);
         m_muon_pair_muon1_mspt        .push_back(m_muon_mspt        [i]);
         m_muon_pair_muon1_msp         .push_back(m_muon_msp         [i]);
         m_muon_pair_muon1_elosstype   .push_back(m_muon_elosstype   [i]);
         m_muon_pair_muon1_trk_pt      .push_back(m_muon_trk_pt      [i]);
         m_muon_pair_muon1_trk_eta     .push_back(m_muon_trk_eta     [i]);
         m_muon_pair_muon1_trk_phi     .push_back(m_muon_trk_phi     [i]);

         m_muon_pair_muon2_pt          .push_back(m_muon_pt          [j]);
         m_muon_pair_muon2_eta         .push_back(m_muon_eta         [j]);
         m_muon_pair_muon2_phi         .push_back(m_muon_phi         [j]);
         m_muon_pair_muon2_pt_precorr  .push_back(m_muon_pt_precorr  [j]);
         m_muon_pair_muon2_eta_precorr .push_back(m_muon_eta_precorr [j]);
         m_muon_pair_muon2_phi_precorr .push_back(m_muon_phi_precorr [j]);
         m_muon_pair_muon2_qual        .push_back(m_muon_qual        [j]);
         m_muon_pair_muon2_deltaP_overP.push_back(m_muon_deltaP_overP[j]);
         m_muon_pair_muon2_d0          .push_back(m_muon_d0          [j]);
         m_muon_pair_muon2_d0_err      .push_back(m_muon_d0_err      [j]);
         m_muon_pair_muon2_z0          .push_back(m_muon_z0          [j]);
         m_muon_pair_muon2_trk_index   .push_back(m_muon_trk_index   [j]);
         m_muon_pair_muon2_mspt        .push_back(m_muon_mspt        [j]);
         m_muon_pair_muon2_msp         .push_back(m_muon_msp         [j]);
         m_muon_pair_muon2_elosstype   .push_back(m_muon_elosstype   [j]);
         m_muon_pair_muon2_trk_pt      .push_back(m_muon_trk_pt      [j]);
         m_muon_pair_muon2_trk_eta     .push_back(m_muon_trk_eta     [j]);
         m_muon_pair_muon2_trk_phi     .push_back(m_muon_trk_phi     [j]);

         //Trigger matching
         #if defined(__ATHENA_21p2__) || defined(__ATHENA_24p2__)
         if(m_use_trigger){
           std::vector<const xAOD::IParticle*> myParticles;
           myParticles.push_back(muon_map[i]);
           myParticles.push_back(muon_map[j]);
           for(auto& trigger_name:m_ListOfDiMuonTriggers){
             if(m_matchTool->match(myParticles,trigger_name,0.1)){
               m_trigger_match_map[trigger_name]->push_back(true);
             }
             else{
               m_trigger_match_map[trigger_name]->push_back(false);
             }
           }
         }
         #endif
       }
     }
   }
   return StatusCode::SUCCESS;
}



//Tracks
void TrigRates::InitTracks(TTree *l_OutTree){
  l_OutTree->Branch("trk_numqual",&m_trk_numqual);
  if(m_store_tracks&Track::StoreBasic){
    l_OutTree->Branch("trk_pt"     ,&m_track_pt);
    l_OutTree->Branch("trk_eta"    ,&m_track_eta);
    l_OutTree->Branch("trk_phi"    ,&m_track_phi);
    l_OutTree->Branch("trk_charge" ,&m_track_charge );
    l_OutTree->Branch("trk_qual"   ,&m_track_quality);
    if((m_store_tracks&Track::StoreDetails)){
      l_OutTree->Branch("trk_z0_wrtPV"       ,&m_track_z0_wrtPV       );
      l_OutTree->Branch("trk_d0"             ,&m_track_d0             );
      l_OutTree->Branch("trk_vz"             ,&m_track_vz             );
      l_OutTree->Branch("trk_Ipix_hits"      ,&m_track_Ipix_hits      );
      l_OutTree->Branch("trk_Ipix_expected"  ,&m_track_Ipix_expected  );
      l_OutTree->Branch("trk_NIpix_hits"     ,&m_track_NIpix_hits     );
      l_OutTree->Branch("trk_NIpix_expected" ,&m_track_NIpix_expected );
      l_OutTree->Branch("trk_sct_hits"       ,&m_track_sct_hits       );
      l_OutTree->Branch("trk_pix_hits"       ,&m_track_pix_hits       );
      l_OutTree->Branch("trk_sct_holes"      ,&m_track_sct_holes      );
      l_OutTree->Branch("trk_pix_holes"      ,&m_track_pix_holes      );
      l_OutTree->Branch("trk_sct_dead"       ,&m_track_sct_dead       );
      l_OutTree->Branch("trk_pix_dead"       ,&m_track_pix_dead       );
      l_OutTree->Branch("trk_sct_shared"     ,&m_track_sct_shared     );
      l_OutTree->Branch("trk_pix_shared"     ,&m_track_pix_shared     );
      l_OutTree->Branch("trk_chi2"           ,&m_track_chi2           );
      l_OutTree->Branch("trk_ndof"           ,&m_track_ndof           );
      l_OutTree->Branch("trk_patternRecoInfo",&m_track_patternRecoInfo);
    }
    if((m_store_tracks&Track::StoreTruthLink) && m_store_truth){
      l_OutTree->Branch("trk_truth_index"    ,&m_trk_truth_index    ); 
      l_OutTree->Branch("trk_truth_prob"     ,&m_trk_truth_prob     ); 
      l_OutTree->Branch("trk_truth_barcode"  ,&m_trk_truth_barcode  ); 
      l_OutTree->Branch("trk_truth_IsPrimary",&m_trk_truth_IsPrimary); 
    }
  }
}

StatusCode TrigRates::ProcessTracks(){
   const xAOD::TrackParticleContainer *l_TrackParticleContainer = nullptr;
   if(evtStore()->retrieve(l_TrackParticleContainer,m_trk_container_key).isFailure()){
     ATH_MSG_ERROR("Could not retrieve TrackParticleContainer with key "<<m_trk_container_key.c_str());
     return StatusCode::FAILURE;
   }
   ATH_MSG_DEBUG("Num Tracks="<<l_TrackParticleContainer->size());

   const xAOD::VertexContainer *l_VertexContainer = nullptr;
   if(evtStore()->retrieve(l_VertexContainer,m_vtx_container_key).isFailure()){
      ATH_MSG_ERROR("Could not retrieve VxContainer with key "<<m_vtx_container_key.c_str());
      return StatusCode::FAILURE;
   }
   const auto* priVtx = *(l_VertexContainer->cbegin());
   if (priVtx->vertexType() != xAOD::VxType::PriVtx) {
     ATH_MSG_WARNING( "First vertex is not of type \"Primary Vertex\"." );
     return StatusCode::FAILURE;
   }
   float z_vtx = priVtx->z();

   m_trk_numqual     .clear();m_trk_numqual.assign(8,0);
   if(m_store_tracks&Track::StoreBasic){
     m_track_pt        .clear();
     m_track_eta       .clear();
     m_track_phi       .clear();
     m_track_charge    .clear();
     m_track_quality   .clear();
     m_track_index_temp.clear();
     if(m_store_tracks&Track::StoreDetails){
       m_track_d0             .clear();
       m_track_z0_wrtPV       .clear();
       m_track_vz             .clear();
       m_track_Ipix_hits      .clear();
       m_track_Ipix_expected  .clear();
       m_track_NIpix_hits     .clear();
       m_track_NIpix_expected .clear();
       m_track_sct_hits       .clear();
       m_track_pix_hits       .clear();
       m_track_sct_holes      .clear();
       m_track_pix_holes      .clear();
       m_track_sct_dead       .clear();
       m_track_pix_dead       .clear();
       m_track_sct_shared     .clear();
       m_track_pix_shared     .clear();
       m_track_chi2           .clear();
       m_track_ndof           .clear();
       m_track_patternRecoInfo.clear();
     }
     if(m_store_truth && (m_store_tracks&Track::StoreTruthLink)){
       m_trk_truth_index_temp1.clear();
       m_trk_truth_index      .clear();
       m_trk_truth_prob       .clear();
       m_trk_truth_barcode    .clear();
       m_trk_truth_IsPrimary  .clear();
     }
   }

   int itrk=0; 
   for(const auto* track : *l_TrackParticleContainer){
     float pt    =track->pt    ();
     float eta   =track->eta   ();
     float phi   =track->phi   ();
     float charge=track->charge();
     int quality=GetTrackQuality(track,priVtx);//MyUtils::TrackQuality(track,z_vtx);

     if(pt>400)                                      m_trk_numqual[0]++;
     if((quality&MyUtils::PP_MIN_BIAS)>0  && pt>400) m_trk_numqual[1]++; 
     if((quality&MyUtils::HI_LOOSE   )>0  && pt>400) m_trk_numqual[2]++; //Added April21 2020
     if((quality&MyUtils::HI_TIGHT   )>0  && pt>400) m_trk_numqual[3]++; //Added April21 2020
     //lines below Added Dec 2023
     if(true)                                        m_trk_numqual[4]++;
     if((quality&MyUtils::PP_MIN_BIAS)>0           ) m_trk_numqual[5]++; 
     if((quality&MyUtils::HI_LOOSE   )>0           ) m_trk_numqual[6]++;
     if((quality&MyUtils::HI_TIGHT   )>0           ) m_trk_numqual[7]++;
     
     if(m_store_tracks&Track::StoreBasic){
       m_track_pt     .push_back(pt);
       m_track_eta    .push_back(eta);
       m_track_phi    .push_back(phi);
       m_track_charge .push_back(charge);
       m_track_quality.push_back(quality);
       m_track_index_temp[track]=itrk;
       if(m_store_tracks&Track::StoreDetails){
         m_track_d0             .push_back(track->d0());
         m_track_z0_wrtPV       .push_back(track->z0()+track->vz() - z_vtx);
         m_track_vz             .push_back(track->vz());
         m_track_Ipix_hits      .push_back(track->auxdata<uint8_t>("numberOfInnermostPixelLayerHits"));
         m_track_Ipix_expected  .push_back(track->auxdata<uint8_t>("expectInnermostPixelLayerHit"));
         m_track_NIpix_hits     .push_back(track->auxdata<uint8_t>("numberOfNextToInnermostPixelLayerHits"));
         m_track_NIpix_expected .push_back(track->auxdata<uint8_t>("expectNextToInnermostPixelLayerHit"));
         m_track_sct_hits       .push_back(track->auxdata<uint8_t>("numberOfSCTHits"));
         m_track_pix_hits       .push_back(track->auxdata<uint8_t>("numberOfPixelHits"));
         m_track_sct_holes      .push_back(track->auxdata<uint8_t>("numberOfSCTHoles"));
         m_track_pix_holes      .push_back(track->auxdata<uint8_t>("numberOfPixelHoles"));
         m_track_sct_dead       .push_back(track->auxdata<uint8_t>("numberOfSCTDeadSensors"));
         m_track_pix_dead       .push_back(track->auxdata<uint8_t>("numberOfPixelDeadSensors"));
         m_track_sct_shared     .push_back(track->auxdata<uint8_t>("numberOfSCTSharedHits"));
         m_track_pix_shared     .push_back(track->auxdata<uint8_t>("numberOfPixelSharedHits"));
         m_track_chi2           .push_back(track->auxdata<float  >("chiSquared"));
         m_track_ndof           .push_back(track->auxdata<float  >("numberDoF"));
         m_track_patternRecoInfo.push_back(track->auxdata<unsigned long>("patternRecoInfo"));
       }
       //http://acode-browser2.usatlas.bnl.gov/lxr-rel20/source/atlas/InnerDetector/InDetValidation/InDetPhysValMonitoring/src/InDetPhysValMonitoringTool.cxx?v=release_20_3_0
       if(m_store_truth && (m_store_tracks&Track::StoreTruthLink) ){
         const                 ElementLink<xAOD::TruthParticleContainer> ptruthContainer=
               (track->auxdata<ElementLink<xAOD::TruthParticleContainer>>("truthParticleLink" ));
         if(ptruthContainer.isValid()){
           const xAOD::TruthParticle *associated_truth=*ptruthContainer;
           m_trk_truth_index_temp1   [associated_truth]=itrk; 
           m_trk_truth_index    .push_back(-1);
           m_trk_truth_prob     .push_back(track->auxdata<float>("truthMatchProbability")); 
           m_trk_truth_barcode  .push_back(associated_truth->barcode());
           m_trk_truth_IsPrimary.push_back(IsPrimaryParticle(associated_truth));
         }
         else{
           m_trk_truth_index    .push_back(-1);
           m_trk_truth_prob     .push_back(0); 
           m_trk_truth_barcode  .push_back(-1);
           m_trk_truth_IsPrimary.push_back(false);
         }
       }
       itrk++;
     }
   }
   return StatusCode::SUCCESS;
}


//PixelTracks
void TrigRates::InitPixTracks(TTree *l_OutTree){
  l_OutTree->Branch("pix_numqual",&m_pixel_numqual);
  if(m_store_pix_tracks&Track::StoreBasic){
    l_OutTree->Branch("pix_pt"       ,&m_pixel_pt);
    l_OutTree->Branch("pix_eta"      ,&m_pixel_eta);
    l_OutTree->Branch("pix_phi"      ,&m_pixel_phi);
    l_OutTree->Branch("pix_charge"   ,&m_pixel_charge );
    l_OutTree->Branch("pix_qual"     ,&m_pixel_quality);
    l_OutTree->Branch("pix_z0_wrtPV" ,&m_pixel_z0_wrtPV       );
    l_OutTree->Branch("pix_d0"       ,&m_pixel_d0             );
    l_OutTree->Branch("pix_vz"       ,&m_pixel_vz             );
    if((m_store_pix_tracks&Track::StoreDetails)){
      l_OutTree->Branch("pix_Ipix_hits"      ,&m_pixel_Ipix_hits      );
      l_OutTree->Branch("pix_Ipix_expected"  ,&m_pixel_Ipix_expected  );
      l_OutTree->Branch("pix_NIpix_hits"     ,&m_pixel_NIpix_hits     );
      l_OutTree->Branch("pix_NIpix_expected" ,&m_pixel_NIpix_expected );
      l_OutTree->Branch("pix_sct_hits"       ,&m_pixel_sct_hits       );
      l_OutTree->Branch("pix_pix_hits"       ,&m_pixel_pix_hits       );
      l_OutTree->Branch("pix_sct_holes"      ,&m_pixel_sct_holes      );
      l_OutTree->Branch("pix_pix_holes"      ,&m_pixel_pix_holes      );
      l_OutTree->Branch("pix_sct_dead"       ,&m_pixel_sct_dead       );
      l_OutTree->Branch("pix_pix_dead"       ,&m_pixel_pix_dead       );
      l_OutTree->Branch("pix_sct_shared"     ,&m_pixel_sct_shared     );
      l_OutTree->Branch("pix_pix_shared"     ,&m_pixel_pix_shared     );
      l_OutTree->Branch("pix_chi2"           ,&m_pixel_chi2           );
      l_OutTree->Branch("pix_ndof"           ,&m_pixel_ndof           );
      l_OutTree->Branch("pix_patternRecoInfo",&m_pixel_patternRecoInfo);
    }
  }
}


StatusCode TrigRates::ProcessPixTracks(){
   const xAOD::TrackParticleContainer *l_TrackParticleContainer = nullptr;
   if(evtStore()->retrieve(l_TrackParticleContainer,m_pix_container_key).isFailure()){
     ATH_MSG_ERROR("Could not retrieve TrackParticleContainer with key "<<m_pix_container_key.c_str());
     return StatusCode::FAILURE;
   }
   ATH_MSG_DEBUG("Num Pixel Tracks="<<l_TrackParticleContainer->size());

   const xAOD::VertexContainer *l_VertexContainer = nullptr;
   if(evtStore()->retrieve(l_VertexContainer,m_vtx_container_key).isFailure()){
      ATH_MSG_ERROR("Could not retrieve VxContainer with key "<<m_vtx_container_key.c_str());
      return StatusCode::FAILURE;
   }
   const auto* priVtx = *(l_VertexContainer->cbegin());
   if (priVtx->vertexType() != xAOD::VxType::PriVtx) {
     ATH_MSG_WARNING( "First vertex is not of type \"Primary Vertex\"." );
     return StatusCode::FAILURE;
   }
   float z_vtx = priVtx->z();

   m_pixel_numqual     .clear();m_pixel_numqual.assign(8,0);
   if(m_store_pix_tracks&Track::StoreBasic){
     m_pixel_pt        .clear();
     m_pixel_eta       .clear();
     m_pixel_phi       .clear();
     m_pixel_charge    .clear();
     m_pixel_quality   .clear();
     m_pixel_d0        .clear();
     m_pixel_z0_wrtPV  .clear();
     m_pixel_vz        .clear();
     if(m_store_pix_tracks&Track::StoreDetails){
       m_pixel_Ipix_hits      .clear();
       m_pixel_Ipix_expected  .clear();
       m_pixel_NIpix_hits     .clear();
       m_pixel_NIpix_expected .clear();
       m_pixel_sct_hits       .clear();
       m_pixel_pix_hits       .clear();
       m_pixel_sct_holes      .clear();
       m_pixel_pix_holes      .clear();
       m_pixel_sct_dead       .clear();
       m_pixel_pix_dead       .clear();
       m_pixel_sct_shared     .clear();
       m_pixel_pix_shared     .clear();
       m_pixel_chi2           .clear();
       m_pixel_ndof           .clear();
       m_pixel_patternRecoInfo.clear();
     }
   }

   for(const auto* track : *l_TrackParticleContainer){
     float pt    =track->pt    ();
     float eta   =track->eta   ();
     float phi   =track->phi   ();
     float charge=track->charge();
     int quality=GetTrackQuality(track,priVtx);//MyUtils::TrackQuality(track,z_vtx);

     if(pt>400)                                      m_pixel_numqual[0]++;
     if((quality&MyUtils::PP_MIN_BIAS)>0  && pt>400) m_pixel_numqual[1]++; 
     if((quality&MyUtils::HI_LOOSE   )>0  && pt>400) m_pixel_numqual[2]++;
     if((quality&MyUtils::HI_TIGHT   )>0  && pt>400) m_pixel_numqual[3]++;
     if(true)                                        m_pixel_numqual[4]++;
     if((quality&MyUtils::PP_MIN_BIAS)>0           ) m_pixel_numqual[5]++; 
     if((quality&MyUtils::HI_LOOSE   )>0           ) m_pixel_numqual[6]++;
     if((quality&MyUtils::HI_TIGHT   )>0           ) m_pixel_numqual[7]++;
     
     if(m_store_pix_tracks&Track::StoreBasic){
       m_pixel_pt      .push_back(pt);
       m_pixel_eta     .push_back(eta);
       m_pixel_phi     .push_back(phi);
       m_pixel_charge  .push_back(charge);
       m_pixel_quality .push_back(quality);
       m_pixel_d0      .push_back(track->d0());
       m_pixel_z0_wrtPV.push_back(track->z0()+track->vz() - z_vtx);
       m_pixel_vz      .push_back(track->vz());
       if(m_store_pix_tracks&Track::StoreDetails){
         m_pixel_Ipix_hits      .push_back(track->auxdata<uint8_t>("numberOfInnermostPixelLayerHits"));
         m_pixel_Ipix_expected  .push_back(track->auxdata<uint8_t>("expectInnermostPixelLayerHit"));
         m_pixel_NIpix_hits     .push_back(track->auxdata<uint8_t>("numberOfNextToInnermostPixelLayerHits"));
         m_pixel_NIpix_expected .push_back(track->auxdata<uint8_t>("expectNextToInnermostPixelLayerHit"));
         m_pixel_sct_hits       .push_back(track->auxdata<uint8_t>("numberOfSCTHits"));
         m_pixel_pix_hits       .push_back(track->auxdata<uint8_t>("numberOfPixelHits"));
         m_pixel_sct_holes      .push_back(track->auxdata<uint8_t>("numberOfSCTHoles"));
         m_pixel_pix_holes      .push_back(track->auxdata<uint8_t>("numberOfPixelHoles"));
         m_pixel_sct_dead       .push_back(track->auxdata<uint8_t>("numberOfSCTDeadSensors"));
         m_pixel_pix_dead       .push_back(track->auxdata<uint8_t>("numberOfPixelDeadSensors"));
         m_pixel_sct_shared     .push_back(track->auxdata<uint8_t>("numberOfSCTSharedHits"));
         m_pixel_pix_shared     .push_back(track->auxdata<uint8_t>("numberOfPixelSharedHits"));
         m_pixel_chi2           .push_back(track->auxdata<float  >("chiSquared"));
         m_pixel_ndof           .push_back(track->auxdata<float  >("numberDoF"));
         m_pixel_patternRecoInfo.push_back(track->auxdata<unsigned long>("patternRecoInfo"));
       }
     }
   }
   return StatusCode::SUCCESS;
}


int TrigRates::GetTrackQuality(const xAOD::TrackParticle* track,const xAOD::Vertex* pv){
  int   n_sct_hits      =track->auxdata<unsigned char>("numberOfSCTHits");
  float d0      = track->d0();
  float z0_wrtPV= track->z0()+track->vz()-pv->z();
  float theta   = track->theta();

  //-------------------------------------------------------------------------------------------------
  bool pass_min_bias=false;
  if (m_trkSelTool_MinBias->accept(*track, pv)) pass_min_bias=true;
  //-------------------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------
  bool pass_hi_loose=false;
  if (m_trkSelTool_HILoose->accept(*track, pv)) pass_hi_loose=true;
  //-------------------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------
  bool pass_hi_loose_additional_SCT_hit=true;
  if(!pass_hi_loose || n_sct_hits<7) pass_hi_loose_additional_SCT_hit=false;
  //-------------------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------
  bool pass_hi_loose_tight_d0_z0=true;
  if(!pass_hi_loose || fabs(d0)>1.0 || fabs(z0_wrtPV*sin(theta))>1.0) pass_hi_loose_tight_d0_z0=false;
  //-------------------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------
  bool pass_hi_loose_tighter_d0_z0=true;
  if(!pass_hi_loose || fabs(d0)>0.5 || fabs(z0_wrtPV*sin(theta))>0.5) pass_hi_loose_tighter_d0_z0=false;
  //-------------------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------
  bool pass_hi_tight=false;
  if (m_trkSelTool_HITight->accept(*track, pv)) pass_hi_tight=true;
  //-------------------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------
  bool pass_hi_tight_loose_d0_z0=true;
  if(pass_hi_tight==false){
    #ifndef __ATHENA_24p2__
    const auto& taccept = m_trkSelTool_HITight->getTAccept();
    #else
    const auto& taccept = m_trkSelTool_HITight->accept(*track, pv);
    #endif
    static const auto d0Index = taccept.getCutPosition("D0");
    static const auto z0Index = taccept.getCutPosition("Z0SinTheta");
    static const auto nCuts = taccept.getNCuts();
    auto cutBitset = taccept.getCutResultBitSet();
    cutBitset |= (1 << d0Index) | (1 << z0Index);
    if(cutBitset.count() != nCuts                   ) pass_hi_tight_loose_d0_z0=false;
    if(fabs(d0)>1.5 || fabs(z0_wrtPV*sin(theta))>1.5) pass_hi_tight_loose_d0_z0=false;
  }
  //-------------------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------
  bool pass_hi_tight_tighter_d0_z0=true;
  if(!pass_hi_tight || fabs(d0)>0.5 || fabs(z0_wrtPV*sin(theta))>0.5) pass_hi_tight_tighter_d0_z0=false;
  //-------------------------------------------------------------------------------------------------

  unsigned short    quality =0;
  if(pass_min_bias                   ) quality+=MyUtils::PP_MIN_BIAS;
  if(pass_hi_loose                   ) quality+=MyUtils::HI_LOOSE;
  if(pass_hi_loose_additional_SCT_hit) quality+=MyUtils::HI_LOOSE_7SCT_HITS;
//if(pass_hi_loose_tight_d0_z0       ) quality+=MyUtils::HI_LOOSE_TIGHT_D0_Z0;
//if(pass_hi_loose_tighter_d0_z0     ) quality+=MyUtils::HI_LOOSE_TIGHTER_D0_Z0;
  if(pass_hi_tight_loose_d0_z0       ) quality+=MyUtils::HI_TIGHT_LOOSE_D0_Z0;
  if(pass_hi_tight                   ) quality+=MyUtils::HI_TIGHT;
//if(pass_hi_tight_tighter_d0_z0     ) quality+=MyUtils::HI_TIGHT_TIGHTER_D0_Z0;
  return quality;
} 


void TrigRates::InitL1TE(TTree *l_OutTree){
  l_OutTree->Branch("L1TE"  ,&m_L1TE  ,"L1TE/F");
  l_OutTree->Branch("L1TE24",&m_L1TE24,"L1TE24/F");
}

StatusCode TrigRates::ProcessL1TE(){
  const xAOD::EnergySumRoI *ptrOnL1te=nullptr;
  if(evtStore()->retrieve(ptrOnL1te,m_L1TE_container_key).isFailure()){
     ATH_MSG_ERROR("Could not retrieve EnergySumRoI with key "<<m_L1TE_container_key.c_str());
     return StatusCode::FAILURE;
   }
   m_L1TE   = ptrOnL1te->energyT()/1000;
   m_L1TE24 = ptrOnL1te->energyTRestricted()/1000;
   return StatusCode::SUCCESS;
}


void TrigRates::InitMET(TTree *l_OutTree){
  l_OutTree->Branch("MET_sumet", &m_MET_sumet);
}

StatusCode TrigRates::ProcessMET(){
   const xAOD::MissingETContainer *l_MissingETContainer = nullptr;
   if(evtStore()->retrieve(l_MissingETContainer,m_met_container_key).isFailure()){
     ATH_MSG_ERROR("Could not retrieve MissingETContainer with key "<<m_met_container_key.c_str());
     return StatusCode::FAILURE;
   }
   m_MET_sumet.clear();
   m_MET_sumet.assign(7,-1.0e9);
   for(const auto& MET: *l_MissingETContainer){
      std::string MET_Name =MET->name();
      float       MET_sumet=MET->sumet();
      
      if     (MET_Name=="EMB" ) m_MET_sumet[0]=MET_sumet;
      else if(MET_Name=="EME" ) m_MET_sumet[1]=MET_sumet;
      else if(MET_Name=="FCAL") m_MET_sumet[2]=MET_sumet;
      else if(MET_Name=="HEC" ) m_MET_sumet[3]=MET_sumet; 
      else if(MET_Name=="PEMB") m_MET_sumet[4]=MET_sumet;
      else if(MET_Name=="PEME") m_MET_sumet[5]=MET_sumet;
      else if(MET_Name=="TILE") m_MET_sumet[6]=MET_sumet;
      //else { ATH_MSG_ERROR("Unknown MET index::"<<MET_Name.c_str());exit(0);}
   } 
   return StatusCode::SUCCESS;
}


Float_t m_pass;
std::vector<float> m_EventWeights;
void TrigRates::InitTruth(TTree *l_OutTree){
  if(m_store_truth & Truth::Basic){
    l_OutTree->Branch("EventWeights" ,&m_EventWeights);
    l_OutTree->Branch("truth_pt"     ,&m_truth_pt);
    l_OutTree->Branch("truth_eta"    ,&m_truth_eta);
    l_OutTree->Branch("truth_phi"    ,&m_truth_phi);
    l_OutTree->Branch("truth_charge" ,&m_truth_charge );
    l_OutTree->Branch("truth_id"     ,&m_truth_id     );
    l_OutTree->Branch("truth_barcode",&m_truth_barcode);
    l_OutTree->Branch("truth_qual"   ,&m_truth_quality);
    if(m_store_truth&Truth::StoreParents) {
      l_OutTree->Branch("truth_parents",&m_truth_parents);
      l_OutTree->Branch("truth_status" ,&m_truth_status);
    }
    if(m_store_tracks >= Track::StoreBasic) m_OutTree->Branch("truth_trk_index" ,&m_truth_trk_index );
    if(m_store_single_muon                ) m_OutTree->Branch("truth_muon_index",&m_truth_muon_index);
  }

  if(m_store_truth & Truth::StoreMuSingle){
    l_OutTree->Branch("truth_muon_pt"     ,&m_truth_muon_pt     );
    l_OutTree->Branch("truth_muon_eta"    ,&m_truth_muon_eta    );
    l_OutTree->Branch("truth_muon_phi"    ,&m_truth_muon_phi    );
    l_OutTree->Branch("truth_muon_ch"     ,&m_truth_muon_ch     );
    l_OutTree->Branch("truth_muon_barcode",&m_truth_muon_barcode);
    l_OutTree->Branch("truth_muon_id"     ,&m_truth_muon_id     );
    l_OutTree->Branch("truth_muon_status" ,&m_truth_muon_status );
  }
  if(m_store_truth & Truth::StoreMuPairs){
    l_OutTree->Branch("truth_mupair_pass" ,&m_pass              , "pass/F");
    l_OutTree->Branch("truth_mupair_asym" ,&m_truth_mupair_asym );
    l_OutTree->Branch("truth_mupair_acop" ,&m_truth_mupair_acop );
    l_OutTree->Branch("truth_mupair_kperp",&m_truth_mupair_kperp);
    l_OutTree->Branch("truth_mupair_pt"   ,&m_truth_mupair_pt   );
    l_OutTree->Branch("truth_mupair_y"    ,&m_truth_mupair_y    );
    l_OutTree->Branch("truth_mupair_phi"  ,&m_truth_mupair_phi  );
    l_OutTree->Branch("truth_mupair_m"    ,&m_truth_mupair_m    );

    l_OutTree->Branch("truth_mupair_pt1"    ,&m_truth_mupair_pt1  );
    l_OutTree->Branch("truth_mupair_eta1"   ,&m_truth_mupair_eta1 );
    l_OutTree->Branch("truth_mupair_phi1"   ,&m_truth_mupair_phi1 );
    l_OutTree->Branch("truth_mupair_ch1"    ,&m_truth_mupair_ch1  );
    l_OutTree->Branch("truth_mupair_bar1"   ,&m_truth_mupair_bar1 );
    l_OutTree->Branch("truth_mupair_id1"    ,&m_truth_mupair_id1  );
    l_OutTree->Branch("truth_mupair_status1",&m_truth_mupair_status1);

    l_OutTree->Branch("truth_mupair_pt2"    ,&m_truth_mupair_pt2  );
    l_OutTree->Branch("truth_mupair_eta2"   ,&m_truth_mupair_eta2 );
    l_OutTree->Branch("truth_mupair_phi2"   ,&m_truth_mupair_phi2 );
    l_OutTree->Branch("truth_mupair_ch2"    ,&m_truth_mupair_ch2  );
    l_OutTree->Branch("truth_mupair_bar2"   ,&m_truth_mupair_bar2 );
    l_OutTree->Branch("truth_mupair_id2"    ,&m_truth_mupair_id2  );
    l_OutTree->Branch("truth_mupair_status2",&m_truth_mupair_status2);
  }
}

StatusCode TrigRates::ProcessTruth(){
   m_truth_pt     .clear();
   m_truth_eta    .clear();
   m_truth_phi    .clear();
   m_truth_charge .clear();
   m_truth_id     .clear();
   m_truth_barcode.clear();
   m_truth_quality.clear();
   m_truth_parents.clear();
   m_truth_status .clear();
   m_truth_trk_index .clear();
   m_truth_muon_index.clear();

   //-------------------------------
   const xAOD::EventInfo* l_EventInfo = nullptr;
   if(evtStore()->retrieve(l_EventInfo,m_EventInfo_key).isFailure()){
      ATH_MSG_ERROR(" Could not retrieve EventInfo with key "<<m_EventInfo_key.c_str());
      return StatusCode::FAILURE;
   }
   m_EventWeights.clear();
   for(float weight:l_EventInfo->mcEventWeights()){
    m_EventWeights.push_back(weight);
   }
   //-------------------------------

   const xAOD::TruthParticleContainer *l_TruthParticleContainer;
   if(evtStore()->retrieve(l_TruthParticleContainer,m_truth_container_key).isFailure()){
     ATH_MSG_ERROR("Could not retrieve TruthParticleContainer with key"<<m_truth_container_key);
     return StatusCode::FAILURE;
   }

   int truth_index=0;
   for(auto track_itr=l_TruthParticleContainer->begin();track_itr!=l_TruthParticleContainer->end();track_itr++){
     auto track=(*track_itr);

     //if the condition below is false, we have to store all patricles and not just the stable ones
     if((m_store_truth & Truth::StoreParents)==0){ 
       if(track->status()!=1) continue;
       if(track->pt()<0.0001 || track->pt()<m_min_pT_Truth) continue;
       if(track->barcode()>=200000 || track->barcode()==0) continue;
     //if(fabs(track->charge())<0.1 || fabs(track->eta())>2.5 ) continue;
       if(fabs(track->charge())<0.1) continue;
     }

     float pt    =track->pt    ();
     float eta   =1000;
     if(fabs(pt)>0.00001) eta   =track->eta   ();
     float phi   =track->phi   ();
     float charge=track->charge();
     int   id    =track->pdgId ();
     int quality =1;
     if(track->isStrangeBaryon() ||track->barcode()<=0) quality=-1;

     std::vector<int> parents;
     for(long unsigned int iparent=0;iparent<track->nParents();iparent++){
       if(!track->parent(iparent)){
         std::cout<<"AAAAAAAAAAAA Parent "<<iparent<<" is missing "<<track->nParents()<<std::endl;
         continue;
       }
       parents.push_back(track->parent(iparent)->barcode());
     }

     m_truth_pt      .push_back(pt);
     m_truth_eta     .push_back(eta);
     m_truth_phi     .push_back(phi);
     m_truth_charge  .push_back(charge);
     m_truth_id      .push_back(id);
     m_truth_barcode .push_back(track->barcode());
     m_truth_quality .push_back(quality);
     m_truth_parents .push_back(parents);
     m_truth_status  .push_back( (track->status()!=1)? false:true);

     if(m_store_tracks & Track::StoreBasic){
       int trk_index=-1;
       if(m_trk_truth_index_temp1.find(track)!=m_trk_truth_index_temp1.end()){
         trk_index=m_trk_truth_index_temp1[track];
         if(m_trk_truth_barcode[trk_index]!=track->barcode()){
           ATH_MSG_ERROR("Barcodes Dont Match");
           return StatusCode::FAILURE;
         }
         m_trk_truth_index[trk_index]=truth_index;
       }
       m_truth_trk_index.push_back(trk_index);
     }

     if(m_store_single_muon){
       int muon_index=-1;
       if(m_muon_truth_index_temp1.find(track)!=m_muon_truth_index_temp1.end()){
         muon_index=m_muon_truth_index_temp1[track];
         if(m_muon_truth_barcode[muon_index]!=track->barcode()){
           ATH_MSG_ERROR("Barcodes Dont Match");
           return StatusCode::FAILURE;
         }
         m_muon_truth_index[muon_index]=truth_index;
       }
       m_truth_muon_index.push_back(muon_index);
     }

     truth_index++;
   }


   if(m_store_truth & Truth::StoreMuSingle){
     m_truth_muon_pt      .clear();
     m_truth_muon_eta     .clear();
     m_truth_muon_phi     .clear();
     m_truth_muon_ch      .clear();
     m_truth_muon_barcode .clear();
     m_truth_muon_id      .clear();
     m_truth_muon_status  .clear();
     for(unsigned int index1=0;index1<m_truth_pt.size();index1++){
       int id1_=fabs(m_truth_id[index1]);
       if(id1_!=13 || m_truth_status[index1]!=true) continue;
       // if((id1_!=13 && id1_!=15 && id1_!=11 && id1_!=211) || m_truth_status[index1]!=true) continue;
       m_truth_muon_pt      .push_back(m_truth_pt     [index1]);
       m_truth_muon_eta     .push_back(m_truth_eta    [index1]);
       m_truth_muon_phi     .push_back(m_truth_phi    [index1]);
       m_truth_muon_ch      .push_back(m_truth_charge [index1]);
       m_truth_muon_barcode .push_back(m_truth_barcode[index1]);
       m_truth_muon_id      .push_back(id1_);
       m_truth_muon_status  .push_back(m_truth_status [index1]);
     }
   }

   if(m_store_truth & Truth::StoreMuPairs){
     m_truth_mupair_asym .clear();
     m_truth_mupair_acop .clear();
     m_truth_mupair_kperp.clear();
     m_truth_mupair_pt   .clear();
     m_truth_mupair_y    .clear();
     m_truth_mupair_phi  .clear();
     m_truth_mupair_m    .clear();

     m_truth_mupair_pt1    .clear();
     m_truth_mupair_eta1   .clear();
     m_truth_mupair_phi1   .clear();
     m_truth_mupair_ch1    .clear();
     m_truth_mupair_bar1   .clear();
     m_truth_mupair_id1    .clear();
     m_truth_mupair_status1.clear();

     m_truth_mupair_pt2    .clear();
     m_truth_mupair_eta2   .clear();
     m_truth_mupair_phi2   .clear();
     m_truth_mupair_ch2    .clear();
     m_truth_mupair_bar2   .clear();
     m_truth_mupair_id2    .clear();
     m_truth_mupair_status2.clear();
  
     const double PI=acos(-1.0);
     m_pass=0;
     for(unsigned int index1=0;index1<m_truth_pt.size();index1++){
       int id1_=fabs(m_truth_id[index1]);
       if(id1_!=13 || m_truth_status[index1]!=true) continue;
       // if((id1_!=13 && id1_!=15 && id1_!=11 && id1_!=211) || m_truth_status[index1]!=true) continue;
       float pt1    =m_truth_pt     [index1];
       float eta1   =m_truth_eta    [index1];
       float phi1   =m_truth_phi    [index1];
       float charge1=m_truth_charge [index1];
       float bar1   =m_truth_barcode[index1];
       float id1    =m_truth_id     [index1];
       for(unsigned int index2=index1+1;index2<m_truth_pt.size();index2++){
         int id2_=fabs(m_truth_id[index2]);
         if(id2_!=id1_ || m_truth_status[index2]!=true) continue;
         float pt2    =m_truth_pt     [index2];
         float eta2   =m_truth_eta    [index2];
         float phi2   =m_truth_phi    [index2];
         float charge2=m_truth_charge [index2];
         float bar2   =m_truth_barcode[index2];
         float id2    =m_truth_id     [index2];
  
         if(id1_==13 && pt1>=3800 && pt2>=3800 && fabs(eta1)<=2.4 && fabs(eta2)<=2.4) m_pass=1;
  
         float asym =(pt1-pt2)/(pt1+pt2);
         float temp =(phi1-phi2+PI);
         float acop =atan2(sin(temp),cos(temp))/PI;
         float kperp=acop*PI*(pt1+pt2)/2.0;
  
         float M=0;
         if     (id1_== 13) M= 105.7  ;
         else if(id1_== 15) M=1776.86 ;
         else if(id1_== 11) M=   0.511;
         else if(id1_==211) M= 139.6  ;
         else{
           ATH_MSG_ERROR(" Unknown Particle "<<id1);
           return StatusCode::FAILURE;
         }
         TLorentzVector M1,M2;
         M1.SetPtEtaPhiM(pt1, eta1, phi1, M);
         M2.SetPtEtaPhiM(pt2, eta2, phi2, M);
         TLorentzVector M3=M1+M2;
  
         m_truth_mupair_asym .push_back(asym);
         m_truth_mupair_acop .push_back(acop);
         m_truth_mupair_kperp.push_back(kperp);
  
         m_truth_mupair_pt .push_back(M3.Pt      ());
         m_truth_mupair_y  .push_back(M3.Rapidity());
         m_truth_mupair_phi.push_back(M3.Phi     ());
         m_truth_mupair_m  .push_back(M3.M       ());
  
         m_truth_mupair_pt1    .push_back(pt1 );
         m_truth_mupair_eta1   .push_back(eta1);
         m_truth_mupair_phi1   .push_back(phi1);
         m_truth_mupair_ch1    .push_back(charge1);
         m_truth_mupair_bar1   .push_back(bar1);
         m_truth_mupair_id1    .push_back(id1);
         m_truth_mupair_status1.push_back(m_truth_status [index1]);
  
         m_truth_mupair_pt2    .push_back(pt2 );
         m_truth_mupair_eta2   .push_back(eta2);
         m_truth_mupair_phi2   .push_back(phi2);
         m_truth_mupair_ch2    .push_back(charge2);
         m_truth_mupair_bar2   .push_back(bar2);
         m_truth_mupair_id2    .push_back(id2);
         m_truth_mupair_status2.push_back(m_truth_status [index2]);
       }
     }
   }

   return StatusCode::SUCCESS;
}

void TrigRates::InitTruthVertex(TTree *l_OutTree){
     l_OutTree->Branch("truth_vtx_z"     ,&m_truth_vtx_z     );
     l_OutTree->Branch("truth_vtx_x"     ,&m_truth_vtx_x     );
     l_OutTree->Branch("truth_vtx_y"     ,&m_truth_vtx_y     );
     l_OutTree->Branch("truth_vtx_t"     ,&m_truth_vtx_t     );
}

StatusCode TrigRates::ProcessTruthVertex(){
   const xAOD::TruthVertexContainer *l_TruthVertexContainer = nullptr;
   if(evtStore()->retrieve(l_TruthVertexContainer,m_truth_vtx_container_key).isFailure()){
      ATH_MSG_ERROR("Could not retrieve TruthVxContainer with key "<<m_truth_vtx_container_key.c_str());
      return StatusCode::FAILURE;
   }
   m_truth_vtx_z     .clear();
   m_truth_vtx_x     .clear();
   m_truth_vtx_y     .clear();
   m_truth_vtx_t     .clear();
   for(const auto* vtx : *l_TruthVertexContainer){
      m_truth_vtx_z     .push_back(vtx->z());
      m_truth_vtx_x     .push_back(vtx->x());
      m_truth_vtx_y     .push_back(vtx->y());
      m_truth_vtx_t     .push_back(vtx->t());
   }
   return StatusCode::SUCCESS;
}



void TrigRates::InitMcEvents(TTree *l_OutTree){
  #if defined(__ATHENA_24p2__)
  l_OutTree->Branch("mc_pt"    , &m_mc_pt    );
  l_OutTree->Branch("mc_px"    , &m_mc_px    );
  l_OutTree->Branch("mc_py"    , &m_mc_py    );
  l_OutTree->Branch("mc_pz"    , &m_mc_pz    );
  l_OutTree->Branch("mc_eta"   , &m_mc_eta   );
  l_OutTree->Branch("mc_phi"   , &m_mc_phi   );
  l_OutTree->Branch("mc_pdgid" , &m_mc_pdgid );
  #endif
}

StatusCode TrigRates::ProcessMcEvents(){
  #if defined(__ATHENA_24p2__)
  m_mc_pt   .clear();
  m_mc_px   .clear();
  m_mc_py   .clear();
  m_mc_pz   .clear();
  m_mc_eta  .clear();
  m_mc_phi  .clear();
  m_mc_pdgid.clear();

  const McEventCollection *mcColl = nullptr;
  if(evtStore()->retrieve(mcColl,"TruthEvent").isFailure()){
     ATH_MSG_ERROR("Could not retrieve McEventCollection with key "<<"TruthEvent");
     return StatusCode::FAILURE;
  }
  ATH_MSG_DEBUG("Retrieved McEventCollection with key "<<"TruthEvent"<<" "<<mcColl->size());


  for (unsigned int cntr = 0; cntr < mcColl->size(); ++cntr) {
    const HepMC::GenEvent* event = (*mcColl)[cntr];
    // Print the event particle/vtx contents
    ATH_MSG_DEBUG("Printing signal event "<<cntr<<" of "<<mcColl->size());

    for (const auto& vertex: event->vertices()){
      for (const auto& particle: vertex->particles_out()) { 
        auto momentum = particle->momentum();
        int  pid      = particle->pdg_id();
        int  status   = particle->status();
        ATH_MSG_DEBUG("Printing particle: "<<status<<" "<<pid<<" "<<momentum.perp()<<" "<<momentum.pseudoRapidity()<<" "<<momentum.phi()<<" "<<particle->generated_mass());
        if(status==1){
          m_mc_pt   .push_back(momentum.perp()); 
          m_mc_px   .push_back(momentum.x()); 
          m_mc_py   .push_back(momentum.y()); 
          m_mc_pz   .push_back(momentum.z()); 
          m_mc_eta  .push_back(momentum.pseudoRapidity()); 
          m_mc_phi  .push_back(momentum.phi()); 
          m_mc_pdgid.push_back(pid); 
        }
      }
    }
  }
  return StatusCode::SUCCESS;
  #else
  ATH_MSG_ERROR("McEventCollection is not implemented for this release");
  return StatusCode::FAILURE;
  #endif
}
