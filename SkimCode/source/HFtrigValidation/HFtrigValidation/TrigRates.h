#ifndef __TrigRates_H__
#define __TrigRates_H__
#include "HFtrigValidation/AthenaVersion.h"

#include "AthenaBaseComps/AthAlgorithm.h"

#include "string"
#include "map"
#include "vector"

#include "GaudiKernel/ToolHandle.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "InDetTrackSelectionTool/IInDetTrackSelectionTool.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODMuon/MuonContainer.h"

#include "PMGAnalysisInterfaces/IPMGTruthWeightTool.h"


// #include "TrigMuonMatching/ITrigMuonMatching.h"
#include "HFtrigValidation/Module_TrigMuonMatching.h"

#if defined(HF_IS_R21) || defined(HF_IS_R25)
  #include "AsgTools/AnaToolHandle.h"
  #include "AsgAnalysisInterfaces/IGoodRunsListSelectionTool.h"
  #include "MuonAnalysisInterfaces/IMuonSelectionTool.h"
  #include "TriggerMatchingTool/IMatchingTool.h"
  #include "MuonAnalysisInterfaces/IMuonCalibrationAndSmearingTool.h"
  #include "MuonAnalysisInterfaces/IMuonEfficiencyScaleFactors.h"
#else
  #include "GoodRunsLists/IGoodRunsListSelectionTool.h"
  #include "MuonSelectorTools/IMuonSelectionTool.h"
  #include "MuonMomentumCorrections/IMuonCalibrationAndSmearingTool.h"
#endif

#if defined(HF_IS_R25)
  #include "ZdcAnalysis/IZdcAnalysisTool.h"
#endif


class TTree;
class Module;

class TrigRates : public AthAlgorithm{
public:
   /** Standard Athena-Algorithm Constructor */
   TrigRates(const std::string& name, ISvcLocator* pSvcLocator);
   /** Default Destructor */
   ~TrigRates() {};


   virtual StatusCode initialize();
   virtual StatusCode execute();
   virtual StatusCode finalize();


private:
   StatusCode CleaningCuts    (bool &pass_cleaning_cuts);


   bool                        m_is_evgen       =false;//=true for EvGen only (withour reconstruction)
   bool                        m_use_trigger    = true;//=true then use the Trigger (make branches, use decision)
   bool                        m_store_MC_weight_names;//=true then store weight names
   bool                        m_store_L1       =false;//set to true for storing L1 TBP/TAP/TAV values
   std::string                 m_Trigger_Chains       ;//Trigger Chains for TrigDecision and branches
   std::string                 m_Muon_Trigger_Chains  ;//Trigger Chains for TrigDecision and branches
   std::string                 m_DiMuon_Trigger_Chains;//Trigger Chains for TrigDecision and branches
   std::vector<std::string>    m_ListOfTriggers       ;
   std::vector<std::string>    m_ListOfMuonTriggers   ;
   std::vector<std::string>    m_ListOfDiMuonTriggers ;

   double                      m_muonTriggerMatchDR;

   std::map<std::string,Bool_t*             > m_trigger_map            ;//flags that indicate if the trigger passes,
   std::map<std::string,Float_t*            > m_trigger_prescale_map   ;//trigger prescale values, 
   std::map<std::string,Bool_t*             > m_trigger_isPrescaled_map;//if a trigger was prescaled in the present event, 
   std::map<std::string,Bool_t*             > m_trigger_L1TBP_map      ;//decision at L1 before prescale, 
   std::map<std::string,Bool_t*             > m_trigger_L1TAP_map      ;//decision at L1 after  prescale, 
   std::map<std::string,Bool_t*             > m_trigger_L1TAV_map      ;//decision at L1 after  veto, 
   std::map<std::string,Bool_t*             > m_trigger_isRerun_map    ;//if Trigger was rerun,
   std::map<std::string,Bool_t*             > m_trigger_rerun_map      ;//Trigger rerun decision,
   
   std::map<std::string,std::vector<Bool_t>*> m_trigger_match_map      ;//matching of reco objects to trigger,
   std::map<std::string,std::vector<Bool_t>*> m_trigger_match_map_0_01 ;//matching of reco objects to trigger, using dR of 0.01
   std::map<std::string,std::vector<Bool_t>*> m_trigger_match_map_V2   ;//matching of reco objects to trigger; doesnot require trigger passing "Physics"
   std::map<std::string,std::vector<Bool_t>*> m_trigger_match_map_V3   ;//matching of reco objects to trigger; require trigger passing "Physics|allowResurrectedDecision"
   
   // Per-leg dimuon matching (split from pair<Bool_t,Bool_t> to avoid ROOT dict issues)
   // leg1 = result.first  (muX matched leg1 in assignment muX->leg1)
   // leg2 = result.second (muX matched leg2 in assignment muX->leg2)
   std::map<std::string,std::vector<Bool_t>*> m_dimu_trigger_match_map_mu1_leg1_0_01;
   std::map<std::string,std::vector<Bool_t>*> m_dimu_trigger_match_map_mu1_leg2_0_01;
   std::map<std::string,std::vector<Bool_t>*> m_dimu_trigger_match_map_mu2_leg1_0_01;
   std::map<std::string,std::vector<Bool_t>*> m_dimu_trigger_match_map_mu2_leg2_0_01;
   std::map<std::string,std::vector<Bool_t>*> m_dimu_trigger_match_map_mu1_leg1_dR;
   std::map<std::string,std::vector<Bool_t>*> m_dimu_trigger_match_map_mu1_leg2_dR;
   std::map<std::string,std::vector<Bool_t>*> m_dimu_trigger_match_map_mu2_leg1_dR;
   std::map<std::string,std::vector<Bool_t>*> m_dimu_trigger_match_map_mu2_leg2_dR;
   
   Bool_t  m_trigger_Flag                     [1000];//stores if a trigger passes, this is written to the output tree
   Float_t m_trigger_prescale_Value           [1000];//stores trigger prescale
   Bool_t  m_trigger_isPrescaled_Flag         [1000];//stores=1 if trigger didnot fire because it was Prescaled
   Bool_t  m_trigger_L1TBP_Flag               [1000];//stores Trigger decision at L1 before prescale
   Bool_t  m_trigger_L1TAP_Flag               [1000];//stores Trigger decision at L1 after  prescale
   Bool_t  m_trigger_L1TAV_Flag               [1000];//stores Trigger decision at L1 after  veto
   Bool_t  m_trigger_isRerun_Flag             [1000];//stores if a trigger was rerun
   Bool_t  m_trigger_rerun_Flag               [1000];//stores if a trigger passes rerun
   
   std::vector<Bool_t> m_trigger_match_Flag     [1000];//whether a given object fired this trigger
   std::vector<Bool_t> m_trigger_match_Flag_0_01[1000];//whether a given object fired this trigger
   std::vector<Bool_t> m_trigger_match_Flag_V2  [1000];//whether a given object fired this trigger; doesnot require trigger passing "Physics" 
   std::vector<Bool_t> m_trigger_match_Flag_V3  [1000];//whether a given object fired this trigger; require trigger passing "Physics|allowResurrectedDecision"
   
   // Per-leg storage: leg1=result.first, leg2=result.second for each of mu1/mu2 at each dR
   // (symmetric chain e.g. HLT_2mu4: leg1==leg2 always; asymmetric e.g. HLT_mu4_mu4noL1: leg1=seeded, leg2=fullscan)
   std::vector<Bool_t> m_dimu_trigger_match_Flag_mu1_leg1_0_01[1000];
   std::vector<Bool_t> m_dimu_trigger_match_Flag_mu1_leg2_0_01[1000];
   std::vector<Bool_t> m_dimu_trigger_match_Flag_mu2_leg1_0_01[1000];
   std::vector<Bool_t> m_dimu_trigger_match_Flag_mu2_leg2_0_01[1000];
   std::vector<Bool_t> m_dimu_trigger_match_Flag_mu1_leg1_dR[1000];
   std::vector<Bool_t> m_dimu_trigger_match_Flag_mu1_leg2_dR[1000];
   std::vector<Bool_t> m_dimu_trigger_match_Flag_mu2_leg1_dR[1000];
   std::vector<Bool_t> m_dimu_trigger_match_Flag_mu2_leg2_dR[1000];
   void       AddTriggerBranches();
   StatusCode ProcessTriggers (bool &pass_trigger_cuts);



   int m_store_EventInfo   =1;//>=1 then store even info
   std::string m_EventInfo_key    ;
   void       InitEventInfo   (TTree *l_OutTree = nullptr);
   StatusCode ProcessEventInfo();
   UInt_t RunNumber=0,lumi_block=0,bunch_crossing_id=0;
   ULong64_t eventNumber=0;
   Float_t AvgIntPerXing=0,ActIntPerXing=0;
   Float_t ZdcEtA=0,ZdcEtC=0;
   Int_t   m_NumTrackNoCuts=0, m_NumTrackPPMinBias=0, m_NumTrackHILoose=0, m_NumTrackHITight=0;
   Float_t m_FCalET,m_FCalETP,m_FCalETN;

   UInt_t m_year; // needed for PbPb centrality

   bool        m_is_Zdc_Calib  = false; //=true only for Zdc Calib Stream
   int         m_store_Zdc     = 0;     //bitflag: 1=basic ZDC, 2=RPD/centroid info, 3=both
   std::string m_ZdcAuxSuffix;          //AuxSuffix for Zdc (e.g. "RP" when reprocessing)
   void       InitZdc   (TTree *l_OutTree);
   StatusCode ProcessZdc();

   bool        m_store_Vtx           = true; //=true then store vertex info           
   std::string m_vtx_container_key;
   void       InitVertex   (TTree *l_OutTree = nullptr);
   StatusCode ProcessVertex();
   std::vector<float> m_vtx_z     ;
   std::vector<float> m_vtx_x     ;
   std::vector<float> m_vtx_y     ;
   std::vector<int>   m_vtx_ntrk  ;


   int         m_store_tracks        =1+2; //1 then store tracks, 2 then store more details, 4 store link to truth
   std::string m_trk_container_key;
   void       InitTracks   (TTree *l_OutTree = nullptr);
   StatusCode ProcessTracks();
   int GetTrackQuality(const xAOD::TrackParticle* track,const xAOD::Vertex* pv);
   std::vector<int>   m_trk_numqual;
   std::vector<float> m_track_pt     ;
   std::vector<float> m_track_eta    ;
   std::vector<float> m_track_phi    ;
   std::vector<float> m_track_charge ;
   std::vector<int>   m_track_quality;
   std::map<const xAOD::TrackParticle*,int> m_track_index_temp;//to link muons to tracks
   std::vector<float> m_track_d0     ;
   std::vector<float> m_track_z0_wrtPV      ;
   std::vector<float> m_track_vz            ;
   std::vector<int>   m_track_Ipix_hits     ;
   std::vector<int>   m_track_Ipix_expected ;
   std::vector<int>   m_track_NIpix_hits    ;
   std::vector<int>   m_track_NIpix_expected;
   std::vector<int>   m_track_sct_hits      ;
   std::vector<int>   m_track_pix_hits      ;
   std::vector<int>   m_track_sct_holes     ;
   std::vector<int>   m_track_pix_holes     ;
   std::vector<int>   m_track_sct_dead      ;
   std::vector<int>   m_track_pix_dead      ;
   std::vector<int>   m_track_sct_shared    ;
   std::vector<int>   m_track_pix_shared    ;
   std::vector<float> m_track_chi2          ;
   std::vector<float> m_track_ndof          ;
   std::vector<unsigned long> m_track_patternRecoInfo;
   std::map<const xAOD::TruthParticle*,int> m_trk_truth_index_temp1;
   std::vector<int  > m_trk_truth_index;
   std::vector<float> m_trk_truth_prob ;
   std::vector<int  > m_trk_truth_barcode;
   std::vector<Bool_t>m_trk_truth_IsPrimary;
   bool IsPrimaryParticle(const xAOD::TruthParticle* particle);
   std::string DoubleToString(double value);

   //-------------------------------------------------------------------
   int         m_store_pix_tracks  =1+2;
   std::string m_pix_container_key;
   void       InitPixTracks   (TTree *l_OutTree = nullptr);
   StatusCode ProcessPixTracks();
   std::vector<int>   m_pixel_numqual       ;
   std::vector<float> m_pixel_pt            ;
   std::vector<float> m_pixel_eta           ;
   std::vector<float> m_pixel_phi           ;
   std::vector<float> m_pixel_charge        ;
   std::vector<int>   m_pixel_quality       ;
   std::vector<float> m_pixel_d0            ;
   std::vector<float> m_pixel_z0_wrtPV      ;
   std::vector<float> m_pixel_vz            ;
   std::vector<int>   m_pixel_Ipix_hits     ;
   std::vector<int>   m_pixel_Ipix_expected ;
   std::vector<int>   m_pixel_NIpix_hits    ;
   std::vector<int>   m_pixel_NIpix_expected;
   std::vector<int>   m_pixel_sct_hits      ;
   std::vector<int>   m_pixel_pix_hits      ;
   std::vector<int>   m_pixel_sct_holes     ;
   std::vector<int>   m_pixel_pix_holes     ;
   std::vector<int>   m_pixel_sct_dead      ;
   std::vector<int>   m_pixel_pix_dead      ;
   std::vector<int>   m_pixel_sct_shared    ;
   std::vector<int>   m_pixel_pix_shared    ;
   std::vector<float> m_pixel_chi2          ;
   std::vector<float> m_pixel_ndof          ;
   std::vector<unsigned long> m_pixel_patternRecoInfo;
   //-------------------------------------------------------------------

   std::string m_met_container_key;
   void       InitMET   (TTree *l_OutTree = nullptr);
   StatusCode ProcessMET();
   std::vector<float> m_MET_sumet;


   bool        m_store_truth_Vtx     = 0; //=true then store truth vertex info           
   std::string m_truth_vtx_container_key;
   void       InitTruthVertex(TTree *l_OutTree = nullptr);
   StatusCode ProcessTruthVertex();
   std::vector<float> m_truth_vtx_z;
   std::vector<float> m_truth_vtx_x;
   std::vector<float> m_truth_vtx_y;
   std::vector<float> m_truth_vtx_t;

   std::string m_truth_event_container_key;

   int         m_store_truth;        //=true then store truth track info
   double      m_min_pT_Truth;
   std::string m_truth_container_key;
   void InitTruth(TTree *l_OutTree = nullptr);
   StatusCode ProcessTruth();
   std::vector<float> m_truth_pt     ;
   std::vector<float> m_truth_eta    ;
   std::vector<float> m_truth_phi    ;
   std::vector<float> m_truth_m    ;
   std::vector<float> m_truth_e    ;
   std::vector<float> m_truth_charge ;
   std::vector<int>   m_truth_id     ;
   std::vector<int>   m_truth_barcode;
   std::vector<int>   m_truth_quality;
   std::vector<int>  m_truth_status ;
   std::vector<std::vector<int>>   m_truth_parents;
   std::vector<std::vector<int>>   m_truth_children;
   std::vector<int>   m_truth_trk_index ;
   std::vector<int>   m_truth_muon_index;

   int m_store_McEvent;
   std::vector<float> m_mc_pt   ;
   std::vector<float> m_mc_px   ;
   std::vector<float> m_mc_py   ;
   std::vector<float> m_mc_pz   ;
   std::vector<float> m_mc_eta  ;
   std::vector<float> m_mc_phi  ;
   std::vector<int  > m_mc_pdgid;
   void InitMcEvents(TTree *l_OutTree);
   StatusCode ProcessMcEvents();

  #if defined(HF_IS_R25)
   std::vector<bool > m_muon_match_mu4roi;
   std::vector<bool > m_muon_match_mu6roi;
  #endif
   
   std::vector<float> m_truth_muon_pt  ;
   std::vector<float> m_truth_muon_eta ;
   std::vector<float> m_truth_muon_phi ;
   std::vector<int> m_truth_muon_ch  ;
   std::vector<int> m_truth_muon_barcode;
   std::vector<int> m_truth_muon_id  ;
   std::vector<int > m_truth_muon_status;

   std::vector<float> m_truth_mupair_asym ;
   std::vector<float> m_truth_mupair_acop ;
   std::vector<float> m_truth_mupair_kperp;
   std::vector<float> m_truth_mupair_pt   ;
   std::vector<float> m_truth_mupair_y    ;
   std::vector<float> m_truth_mupair_phi  ;
   std::vector<float> m_truth_mupair_m    ;

   std::vector<float> m_truth_mupair_pt1  ;
   std::vector<float> m_truth_mupair_eta1 ;
   std::vector<float> m_truth_mupair_phi1 ;
   std::vector<int> m_truth_mupair_ch1  ;
   std::vector<int> m_truth_mupair_bar1 ;
   std::vector<int> m_truth_mupair_id1  ;
   std::vector<int > m_truth_mupair_status1;
   
   std::vector<float> m_truth_mupair_pt2  ;
   std::vector<float> m_truth_mupair_eta2 ;
   std::vector<float> m_truth_mupair_phi2 ;
   std::vector<int> m_truth_mupair_ch2  ;
   std::vector<int> m_truth_mupair_bar2 ;
   std::vector<int> m_truth_mupair_id2  ;
   std::vector<int > m_truth_mupair_status2;


   bool        m_store_L1TE            = true; //=true then store ET info (from MET object)
   std::string m_L1TE_container_key;
   void       InitL1TE   (TTree *l_OutTree = nullptr);
   StatusCode ProcessL1TE();
   Float_t m_L1TE,m_L1TE24;


   std::string m_muons_key;
   std::string m_hlt_muons_key;
   std::string m_hlt_muons_fs_key;
   bool m_store_single_muon;
   bool m_store_acoplanar_muon;
   bool m_store_dimuon_perleg;   // save per-leg branch vectors (default false)
   bool m_store_muon_truth;
   void InitMuons(TTree *l_OutTree = nullptr);
   void ClearMuons();
   StatusCode ProcessMuons();
   std::vector<float> m_muon_pt     ;
   std::vector<float> m_muon_eta    ;
   std::vector<float> m_muon_phi    ;
   std::vector<float> m_muon_pt_precorr ;
   std::vector<float> m_muon_eta_precorr;
   std::vector<float> m_muon_phi_precorr;
   std::vector<int>   m_muon_qual   ;
   std::vector<float> m_muon_deltaP_overP;
   std::vector<float> m_muon_eff_corr_medium_data;
   std::vector<float> m_muon_eff_corr_medium_mc  ;
   std::vector<float> m_muon_eff_SF_medium       ;
   std::vector<float> m_muon_eff_corr_tight_data ;
   std::vector<float> m_muon_eff_corr_tight_mc   ;
   std::vector<float> m_muon_eff_SF_tight        ;
   std::vector<float> m_muon_d0     ;
   std::vector<float> m_muon_d0_err ;
   std::vector<float> m_muon_z0     ;
   std::vector<float> m_muon_trk_p  ;
   std::vector<float> m_muon_trk_pt ;
   std::vector<float> m_muon_trk_eta;
   std::vector<float> m_muon_trk_phi;
   std::vector<int>   m_muon_trk_charge;
   std::vector<int>   m_muon_trk_index;
   std::vector<float> m_muon_trk_vx  ;
   std::vector<float> m_muon_trk_vy  ;
   std::vector<float> m_muon_trk_vz  ;
   std::vector<float> m_muon_me_p   ;
   std::vector<int>   m_muon_elosstype;
   std::vector<float> m_muon_mspt;
   std::vector<float> m_muon_msp ;
   std::vector<int  > m_muon_ms_phi_hits ;
   std::vector<int  > m_muon_ms_eta_hits ;
   std::vector<Bool_t>m_muon_trig_match;
   std::map<const xAOD::TruthParticle*,int> m_muon_truth_index_temp1;
   std::vector<int  > m_muon_truth_index;
   std::vector<float> m_muon_truth_prob ;
   std::vector<int  > m_muon_truth_barcode;
   std::vector<Bool_t>m_muon_truth_IsPrimary;
   std::vector<float> m_muon_truth_pt       ;
   std::vector<float> m_muon_truth_eta      ;
   std::vector<float> m_muon_truth_phi      ;
   std::vector<int  > m_muon_truth_charge   ;
   std::vector<int  > m_muon_truth_id       ;
   std::vector<int  > m_muon_truth_quality  ;
   std::vector<Bool_t>m_muon_truth_status   ; 

   std::vector<float> m_muon_pair_acop         ;
   std::vector<float> m_muon_pair_d0           ;
   std::vector<int>   m_muon_pair_muon1_index  ;
   std::vector<int>   m_muon_pair_muon2_index  ;
   std::vector<float> m_muon_pair_muon1_pt     ;
   std::vector<float> m_muon_pair_muon1_eta    ;
   std::vector<float> m_muon_pair_muon1_phi    ;
   std::vector<float> m_muon_pair_muon1_pt_precorr ;
   std::vector<float> m_muon_pair_muon1_eta_precorr;
   std::vector<float> m_muon_pair_muon1_phi_precorr;
   std::vector<int>   m_muon_pair_muon1_qual   ;
   std::vector<float> m_muon_pair_muon1_deltaP_overP;
   //std::vector<Bool_t>m_muon_pair_muon1_trig_match;
   std::vector<float> m_muon_pair_muon1_d0     ;
   std::vector<float> m_muon_pair_muon1_d0_err ;
   std::vector<float> m_muon_pair_muon1_z0     ;
   std::vector<int>   m_muon_pair_muon1_trk_index;
   std::vector<float> m_muon_pair_muon1_mspt;
   std::vector<float> m_muon_pair_muon1_msp ;
   std::vector<int>   m_muon_pair_muon1_elosstype;
   std::vector<float> m_muon_pair_muon1_trk_pt ;
   std::vector<float> m_muon_pair_muon1_trk_eta;
   std::vector<float> m_muon_pair_muon1_trk_phi;
   std::vector<int>   m_muon_pair_muon1_trk_charge;
   std::vector<float> m_muon_pair_muon2_pt     ;
   std::vector<float> m_muon_pair_muon2_eta    ;
   std::vector<float> m_muon_pair_muon2_phi    ;
   std::vector<float> m_muon_pair_muon2_pt_precorr ;
   std::vector<float> m_muon_pair_muon2_eta_precorr;
   std::vector<float> m_muon_pair_muon2_phi_precorr;
   std::vector<int>   m_muon_pair_muon2_qual   ;
   std::vector<float> m_muon_pair_muon2_deltaP_overP;
   //std::vector<Bool_t>m_muon_pair_muon2_trig_match;
   std::vector<float> m_muon_pair_muon2_d0     ;
   std::vector<float> m_muon_pair_muon2_d0_err ;
   std::vector<float> m_muon_pair_muon2_z0     ;
   std::vector<int>   m_muon_pair_muon2_trk_index;
   std::vector<float> m_muon_pair_muon2_mspt;
   std::vector<float> m_muon_pair_muon2_msp ;
   std::vector<int>   m_muon_pair_muon2_elosstype;
   std::vector<float> m_muon_pair_muon2_trk_pt ;
   std::vector<float> m_muon_pair_muon2_trk_eta;
   std::vector<float> m_muon_pair_muon2_trk_phi;
   std::vector<int>   m_muon_pair_muon2_trk_charge;



   std::vector<Module*> m_modules;
   std::string m_HIEventShapeContainer_key;
   std::vector<std::string> m_track_jet_container_keys;


   //Alg Properties
   bool        m_use_GRL        = true ;     //=true then use GRL
   float       m_MaxZvtx        = 250.f;     //max zvtx in mm
   bool        m_StoreAllEvents = true;     //=true then store all events,=false then store events that pass atleast one of the configured triggers


   //Tree and variables to be written to tree and related stuff
   TTree *m_OutTree = nullptr;


   //Tools
   ToolHandle<IGoodRunsListSelectionTool>      m_grlTool;
   ToolHandle<InDet::IInDetTrackSelectionTool> m_trkSelTool_HITight;
   ToolHandle<InDet::IInDetTrackSelectionTool> m_trkSelTool_HILoose;
   ToolHandle<InDet::IInDetTrackSelectionTool> m_trkSelTool_MinBias;
   ToolHandle<CP::IMuonSelectionTool>          m_muonSelection;
   ToolHandle<CP::IMuonCalibrationAndSmearingTool> m_muonCalibrationAndSmearingTool;
   ToolHandle<CP::IMuonEfficiencyScaleFactors    > m_effi_SF_tool_medium;
   ToolHandle<CP::IMuonEfficiencyScaleFactors    > m_effi_SF_tool_tight ;
   bool m_use_effi_SF_tool = true;
   
   ToolHandle<PMGTools::IPMGTruthWeightTool> m_weightTool;    // weight tool (recording both the name and value of the weights)
   
   #if defined(HF_IS_R25)
     asg::AnaToolHandle<Trig::TrigDecisionTool>  m_trigTool;
     ToolHandle<Trig::IMatchingTool>             m_matchTool;
     // ToolHandle<Trig::ITrigMuonMatching>         m_muonmatchTool; // tool for (di)muon trigger matching
   #elif defined(HF_IS_R21)
     ToolHandle<Trig::TrigDecisionTool>          m_trigTool;
     ToolHandle<Trig::IMatchingTool>             m_matchTool;
     // ToolHandle<Trig::ITrigMuonMatching>         m_muonmatchTool; // tool for (di)muon trigger matching
   #else
     ToolHandle<Trig::ITrigMuonMatching>         m_matchTool;
   #endif
  bool m_ApplyMuonCalibrations=false;
  bool m_isRun3=true;

  // Dimuon trigger matching module (owned by this algorithm)
  TrigMuonMatchingModule* m_dimuMatchModule = nullptr;

};
#endif
