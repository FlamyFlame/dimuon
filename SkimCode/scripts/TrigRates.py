
#data with muons
do_pp2015 = False
do_pp2017 = False
do_pp2024 = False
do_hi2018 = False
do_hi2015 = False
do_hi2023 = False
do_hi2024 = True

#overlay
do_pp_MC_fullsim_17       =False
do_pp_MC_fullsim_24       =False



#------------------------------------------------------------
#wheather the sample is HION or not
is_HION=False
if(do_hi2018     or 
   do_hi2015     or 
   do_hi2023     or
   do_hi2024
   ):
  is_HION=True
#------------------------------------------------------------


#------------------------------------------------------------
#for Run3 data
is_Run3=False
if (do_hi2023 or do_hi2024 or do_pp2024 or do_pp_MC_fullsim_24):
  is_Run3=True
#------------------------------------------------------------


#------------------------------------------------------------
#for Analysis in Releases>=22 and (athena) Releases>=24
is_Rel22orAbove=False
is_Rel24orAbove=False
if (do_hi2023 or do_hi2024 or do_pp2024 or do_pp_MC_fullsim_24):
  is_Rel22orAbove=True
if (do_hi2024 or do_pp2024 or do_pp_MC_fullsim_24):
  is_Rel24orAbove=True

#------------------------------------------------------------

#------------------------------------------------------------
is_MC = do_pp_MC_fullsim_17 or do_pp_MC_fullsim_24
#------------------------------------------------------------

#------------------------------------------------------------
RunYear         =0
m_EvtMax        =1000
dataSource      ='data'
GRL             =[]
MinBias_triggers=[]
Muon_triggers   =[]
DiMuon_triggers =[]

if do_hi2018:
  RunYear          =2018 # only set RunYear for PbPb (needed for centrality)
  GRL              =["data18_hi.periodAllYear_DetStatus-v106-pro22-14_Unknown_PHYS_HeavyIonP_All_Good.xml"]
  InputFile        ="/afs/cern.ch/user/s/soumya/workarea/DATA/JobTestData/data18_hi/data18_hi.00366142.physics_HardProbes.merge.AOD.f1027_m2037._lb0570._0003.1"
  MinBias_triggers =[]
  DiMuon_triggers  =["HLT_2mu3","HLT_2mu4","HLT_2mu6","HLT_2mu8","HLT_mu4_mu4noL1"]
elif do_hi2015:
  RunYear          =2015 # only set RunYear for PbPb (needed for centrality)
  GRL              =["data15_hi.periodAllYear_DetStatus-v75-repro20-01_DQDefects-00-02-02_PHYS_HeavyIonP_All_Good.xml"]
  InputFile        ="/afs/cern.ch/user/s/soumya/workarea/DATA/JobTestData/data15_hi/AOD.16615737._000065.pool.root.1"
  MinBias_triggers =[]
  DiMuon_triggers  =["HLT_2mu3","HLT_2mu4","HLT_2mu6","HLT_2mu8","HLT_mu4_mu4noL1"]
elif do_hi2023:
  RunYear          =2023 # only set RunYear for PbPb (needed for centrality)
  GRL              =["data23_hi.periodAllYear_DetStatus-v113-pro31-08_MERGED_PHYS_HeavyIonP_All_Good.xml"]
  InputFile        ="/eos/user/y/yuhang/data/data23_hi_testfile_AOD/data23_hi.00462240.physics_HardProbes.merge.AOD.f1399_m2209._lb0409._0002.1"
  Muon_triggers    =["HLT_mu4_L1MU3V"  ,
                     "HLT_mu6_L1MU3V"  ,
                     "HLT_mu6_L1MU5VF" ,
                     "HLT_mu4_L1MU3V_VTE50",
                     "HLT_mu6_L1MU3V_VTE50",
                     "HLT_mu8_L1MU5VF_VTE50"]
  DiMuon_triggers  =["HLT_2mu4_L12MU3V","HLT_mu4_mu4noL1_L1MU3V"]
elif do_hi2024:
  RunYear          =2024 # only set RunYear for PbPb (needed for centrality)
  GRL              =["physics_HI2024_50ns.xml"]
  InputFile        ="/eos/user/y/yuhang/data/data24_hi_testfile_AOD/data24_hi.00489961.physics_HardProbes.merge.AOD.f1550_m2267._lb0113._0007.1"
  Muon_triggers    =["HLT_mu4_L1MU3V"  ,
                     "HLT_mu6_L1MU3V"  ,
                     "HLT_mu6_L1MU5VF" ,
                     "HLT_mu8_L1MU5VF" ,
                     "HLT_mu10_L1MU8F" ,
                     "HLT_mu10_L1MU5VF",
                     "HLT_mu4noL1_hi_uccTh3_L1jTE10000", 
                     "HLT_mu4noL1_hi_uccTh2_L1jTE9000" , 
                     "HLT_mu4noL1_L1ZDC_HELT25_jTE4000", 
                     "HLT_mu4noL1_L1ZDC_HELT20_jTE4000", 
                     "HLT_mu4noL1_L1ZDC_HELT15_jTE4000"]
  DiMuon_triggers  =["HLT_2mu4_L12MU3V", "HLT_mu4_mu4noL1_L1MU3V"]
elif do_pp2024:
  GRL              =["physics_2024ppRef_25ns.xml"]
  InputFile        ="/eos/user/y/yuhang/data/data24_pp_testfile_AOD/data24_5p36TeV.00488427.physics_Main.merge.AOD.f1529_m2259._lb0128._0004.1"
  Muon_triggers    =["HLT_mu4_L1MU3V"  ,
                     "HLT_mu6_L1MU3V"  ,
                     "HLT_mu8_L1MU5VF" ,
                     "HLT_mu10_L1MU8F",
                     "HLT_mu12_L1MU8F",
                     "HLT_mu15_L1MU8F",
                     "HLT_mu15_L1MU14FCH"]
  DiMuon_triggers  =["HLT_2mu3_L12MU3V","HLT_2mu4_L12MU3V","HLT_mu4_mu6_L12MU3V","HLT_mu4_mu4noL1_L1MU3V"]
elif do_pp2017:
  GRL              =["data17_5TeV.periodAllYear_DetStatus-v98-pro21-16_Unknown_PHYS_StandardGRL_All_Good_25ns_ignore_GLOBAL_LOWMU.xml"]
  InputFile        ="/eos/user/y/yuhang/data/pp_17/data17_5TeV/*"
  Muon_triggers    =["HLT_mu4", "HLT_mu6", "HLT_mu8", "HLT_mu10"]
  DiMuon_triggers  =["HLT_2mu4"]
elif do_pp2015:
  GRL              =["data15_5TeV.periodAllYear_DetStatus-v75-repro20-01_DQDefects-00-02-02_PHYS_HeavyIonP_All_Good.xml",
                     "data15_5TeV.periodVdM_DetStatus-v75-repro20-01_DQDefects-00-02-02_PHYS_HeavyIonP_All_Good.xml"]
  # InputFile        ="/eos/user/y/yuhang/data/pp_15/data15_5TeV/*"
  InputFile        ="/eos/user/y/yuhang/data/pp_15/data15_5TeV_r9582/*"
  Muon_triggers    =["HLT_mu4", "HLT_mu6", "HLT_mu8", "HLT_mu10"]
  DiMuon_triggers  =["HLT_2mu4"]
elif do_pp_MC_fullsim_17:
  dataSource       ='geant4'
  InputFile        ="/eos/user/y/yuhang/data/mc_powheg_fullsim/example_fullsim_AOD_files/AOD.36949837._016835.pool.root.1"
else :
  print("*"*50,"\nUnknown Dataset\n","*"*50)
  exit()
print(MinBias_triggers+Muon_triggers+DiMuon_triggers)
#------------------------------------------------------------








import AthenaPoolCnvSvc.ReadAthenaPool
import glob
svcMgr.EventSelector.InputCollections = glob.glob(InputFile)
theApp.EvtMax = m_EvtMax


from AthenaCommon.GlobalFlags import globalflags
globalflags.DataSource.set_Value_and_Lock(dataSource) #data geant4
svcMgr += CfgMgr.AthenaEventLoopMgr(EventPrintoutInterval=1)
globalflags.DatabaseInstance = 'CONDBR2' # e.g. if you want run 2 data,  set to COMP200 if you want run1. This is completely ignored in MC.
#globalflags.DatabaseInstance = 'COMP200' # e.g. if you want run 2 data,  set to COMP200 if you want run1. This is completely ignored in MC.
#TODO:: WHAT ABOUT RUN3??


#the tool service
from AthenaCommon.AppMgr import ToolSvc

##------------------------------------------------------------
#Configure Track selection tool
#https://twiki.cern.ch/twiki/bin/view/AtlasProtected/InDetTrackSelectionTool
ToolSvc += CfgMgr.InDet__InDetTrackSelectionTool("TrackSelectionTool_MinBias",CutLevel="MinBias",minPt=100.)
ToolSvc += CfgMgr.InDet__InDetTrackSelectionTool("TrackSelectionTool_HILoose",CutLevel="HILoose",minPt=100.)
ToolSvc += CfgMgr.InDet__InDetTrackSelectionTool("TrackSelectionTool_HITight",CutLevel="HITight",minPt=100.)
##------------------------------------------------------------

##---------------------------- Muon Selection Tool --------------------------------
# Create a MuonSelectionTool if we do not yet have one 
if not hasattr(ToolSvc, "MyMuonSelectionTool"):
    from MuonSelectorTools.MuonSelectorToolsConf import CP__MuonSelectionTool
    ToolSvc += CP__MuonSelectionTool("MyMuonSelectionTool")
    ToolSvc.MyMuonSelectionTool.MaxEta = 2.5

    if(is_HION):
      ToolSvc.MyMuonSelectionTool.TrtCutOff=True

    if is_Rel22orAbove:
      ToolSvc.MyMuonSelectionTool.IsRun3Geo=False
      if(is_Run3):
        ToolSvc.MyMuonSelectionTool.IsRun3Geo=True
    print(ToolSvc.MyMuonSelectionTool)
##------------------------------------------------------------

##---------------------------- Muon Calibration Tool --------------------------------
if is_Rel22orAbove:
  #https://gitlab.cern.ch/atlas/athena/-/blob/release/24.2.9/PhysicsAnalysis/MuonID/MuonIDAnalysis/MuonMomentumCorrections/MuonMomentumCorrections/MuonCalibTool.h
  #https://twiki.cern.ch/twiki/bin/view/AtlasProtected/MCPAnalysisGuidelinesR22#CP_MuonCalibTool_tool
  ToolSvc += CfgMgr.CP__MuonCalibTool("MyMuonCalibrationAndSmearingTool")
  ToolSvc.MyMuonCalibrationAndSmearingTool.calibMode=0
  ToolSvc.MyMuonCalibrationAndSmearingTool.IsRun3Geo=False
  if(is_Run3):
    ToolSvc.MyMuonCalibrationAndSmearingTool.IsRun3Geo=True
else:
  #https://twiki.cern.ch/twiki/bin/view/AtlasProtected/MCPAnalysisWinterMC16
  #https://twiki.cern.ch/twiki/bin/view/AtlasProtected/MuonMomentumCorrectionsSubgroup
  #At some point the twiki said to use "Year=Data16" for 2015 data, but now it does not
  from MuonMomentumCorrections.MuonMomentumCorrectionsConf import CP__MuonCalibrationAndSmearingTool
  ToolSvc += CP__MuonCalibrationAndSmearingTool("MyMuonCalibrationAndSmearingTool")
  print(ToolSvc.MyMuonCalibrationAndSmearingTool)
  
  ToolSvc.MyMuonCalibrationAndSmearingTool.Year                 ="Data18"
  ToolSvc.MyMuonCalibrationAndSmearingTool.StatComb             =False
  # SagittaRelease is gives the directory in which the correction factors for this release + dataset are stored
  # Regarding sagittaBiasDataAll_03_02_19:
    # R21 corrections: full run 2 with uniform treatment of 2015+2016, 2017, 2018 data.
    # Hence, should not need to change the Data18 substring SagittaRelease
    # it includes MC phase-space correction, 30x30 bins for good stat, iterations set 4, for minimal RMS
  ToolSvc.MyMuonCalibrationAndSmearingTool.SagittaRelease       ="sagittaBiasDataAll_03_02_19_Data18"
  # MuonCalibrationAndSmearingTool in Athena 21.2 and AthAnalysis 21.2 has the property sagittaMapsInputType
  # which has value NOMINAL (0) and SINGLE (1)
  # default value is SINGLE
  # this specifies number of root files (iterations) containing correction factors, which is SagittaRelease-dependent
  # SINGLE means only one root file (one iteration of correction), which is true only for a few specific releases
  # We use the default release sagittaBiasDataAll_03_02_19, which requires NOMINAL input type
  # see: https://acode-browser1.usatlas.bnl.gov/lxr/source/athena/PhysicsAnalysis/MuonID/MuonIDAnalysis/MuonMomentumCorrections/MuonMomentumCorrections/MuonCalibrationAndSmearingTool.h?v=21.2
  # When is_Rel22orAbove is false (Athena 21), we use AthAnalysis 21.2: hence we always need to set the sagittaMapsInputType value to be NOMINAL (0)
  ToolSvc.MyMuonCalibrationAndSmearingTool.sagittaMapsInputType = 0 # MCAST::SagittaInputHistType::NOMINAL
  ToolSvc.MyMuonCalibrationAndSmearingTool.SagittaCorr          =True
  ToolSvc.MyMuonCalibrationAndSmearingTool.doSagittaMCDistortion=False
  ToolSvc.MyMuonCalibrationAndSmearingTool.Release              ="Recs2020_03_03"
  ToolSvc.MyMuonCalibrationAndSmearingTool.SagittaCorrPhaseSpace=True
  
  if   do_hi2015:
    ToolSvc.MyMuonCalibrationAndSmearingTool.Year               ="Data16"
  elif do_hi2018:
    pass
  elif do_pp2017:
    ToolSvc.MyMuonCalibrationAndSmearingTool.Year               ="Data17"
    ToolSvc.MyMuonCalibrationAndSmearingTool.SagittaRelease     ="sagittaBiasDataAll_03_02_19_Data17"
  elif do_pp2015:
    ToolSvc.MyMuonCalibrationAndSmearingTool.Year               ="Data16"
    ToolSvc.MyMuonCalibrationAndSmearingTool.SagittaRelease     ="sagittaBiasDataAll_03_02_19_Data17" # sagittaBiasDataAll_03_02_19_Data15 NOT FOUND
  elif do_pp_MC_fullsim_17:
    ToolSvc.MyMuonCalibrationAndSmearingTool.Year               ="Data17"
    ToolSvc.MyMuonCalibrationAndSmearingTool.SagittaRelease     ="sagittaBiasDataAll_03_02_19_Data17"
  elif do_pp_MC_fullsim_24:
    ToolSvc.MyMuonCalibrationAndSmearingTool.SagittaCorr        =False
  else:
    print("*"*50,"\nCouldnot configure MuonCalibrationAndSmearingTool\n","*"*50)
    exit()
  
print(ToolSvc.MyMuonCalibrationAndSmearingTool)
##------------------------------------------------------------

##---------------------------- Muon Efficiency Tool --------------------------------
#https://gitlab.cern.ch/atlas/athena/blob/21.2/PhysicsAnalysis/MuonID/MuonIDAnalysis/MuonEfficiencyCorrections/Root/MuonEfficiencyScaleFactors.cxx
#https://twiki.cern.ch/twiki/bin/view/Atlas/MuonEfficiencyScaleFactorsToolUsage
if not is_Rel24orAbove:
  from MuonEfficiencyCorrections.CommonToolSetup import *
  scale_medium = GetMuonEfficiencyTool("Medium",Release="210222_Precision_r21") # release R21 recommendations
  scale_tight  = GetMuonEfficiencyTool("Tight" ,Release="210222_Precision_r21")
else:
  scale_medium = None
  scale_tight = None
##------------------------------------------------------------

##---------------------------- Muon Trigger Matching Tool --------------------------------
# Create trigger matching tool
#https://twiki.cern.ch/twiki/bin/view/Atlas/XAODMatchingTool
#https://twiki.cern.ch/twiki/bin/view/Atlas/R22TriggerAnalysis
if not is_MC:
  if is_Run3:
    ToolSvc += CfgMgr.Trig__R3MatchingTool("MyTriggerMatchTool", OutputLevel=DEBUG)
  else:
    ToolSvc += CfgMgr.Trig__MatchingTool("MyTriggerMatchTool", OutputLevel=DEBUG)
##------------------------------------------------------------

##---------------------------- GRL Tool --------------------------------
#The GRLTool
MyGRLTool=CfgMgr.GoodRunsListSelectionTool("MyGRLTool")                                                             
MyGRLTool.GoodRunsListVec=GRL
#MyGRLTool.PassThrough    =True                                                                                      
ToolSvc += MyGRLTool  
##------------------------------------------------------------


if not is_Rel24orAbove:
  MessageSvc.defaultLimit = 1000

#the main algorithm sequence
algSeq = CfgMgr.AthSequencer("AthAlgSeq")


#------------------------------------------------------------
#main algorithm to run
TrigRatesAlg=CfgMgr.TrigRates("TrigRatesAlg",OutputLevel=INFO)
#------------------------------------------------------------
TrigRatesAlg.TrackSelectionTool_MinBias  =ToolSvc.TrackSelectionTool_MinBias
TrigRatesAlg.TrackSelectionTool_HILoose  =ToolSvc.TrackSelectionTool_HILoose
TrigRatesAlg.TrackSelectionTool_HITight  =ToolSvc.TrackSelectionTool_HITight
TrigRatesAlg.MuonSelectionTool           =ToolSvc.MyMuonSelectionTool
if not is_MC:
  TrigRatesAlg.TriggerMatchTool            =ToolSvc.MyTriggerMatchTool

TrigRatesAlg.GRLTool                     =ToolSvc.MyGRLTool
TrigRatesAlg.muonCorrectionTool          =ToolSvc.MyMuonCalibrationAndSmearingTool
TrigRatesAlg.use_effi_SF_tool            = (not is_Rel24orAbove)
TrigRatesAlg.effi_SF_tool_medium         =scale_medium
TrigRatesAlg.effi_SF_tool_tight          =scale_tight
#------------------------------------------------------------

TrigRatesAlg.RunYear                 =RunYear
TrigRatesAlg.IsEvgen                 =False                     #True only for truth-only input (i.e. reco is missing)
TrigRatesAlg.UseGRL                  =True                      #Automatically ignored for MC
TrigRatesAlg.StoreAllEvents          =False                     #Typically true only for MC & CalibStream, ignores Trigger requirement if set to "True"
TrigRatesAlg.UseTrigger              =(not is_MC)               #Typically True for data False for MC; If "False" Trigger branches are not created
TrigRatesAlg.StoreL1Decision         =(not is_MC)               #Typically 0; set to 1 for storing L1 TBP/TAP/TAV values
TrigRatesAlg.TriggerChains           ="|".join(MinBias_triggers)
TrigRatesAlg.MuonTriggerChains       ="|".join(Muon_triggers)
TrigRatesAlg.DiMuonTriggerChains     ="|".join(DiMuon_triggers)
TrigRatesAlg.StoreEventInfo          =1                         
TrigRatesAlg.StoreTracks             =0                         #0+1+2+4+8 bitmap
TrigRatesAlg.StorePixTracks          =0                         #0+1+2+4 bitmap
TrigRatesAlg.MaxZvtx                 =250                       #require vtx with |z_vtx|<MaxZvtx; -negative to disable vtx cut
TrigRatesAlg.StoreVtx                =True                      #Always ON
TrigRatesAlg.StoreL1TE               =(not is_MC)                      #ON only for checks 
TrigRatesAlg.StoreMuonTruth          =is_MC                     #True only for MC
TrigRatesAlg.StoreSingleMuon         =True                      #Typically ON
TrigRatesAlg.StoreAcoplanarMuon      =True                      #OFF for pp Ridge analysis
TrigRatesAlg.HIEventShapeContainerKey="HIEventShape"            #"HIEventShape" or "CaloSums";="" to turn OFF, for exmple in pp
TrigRatesAlg.METContainerKey         =""                        #"MET_Calo";="" to turn OFF
TrigRatesAlg.StoreTruthVtx           =is_MC                     #True only for MC
TrigRatesAlg.StoreTruth              =1+2+4+8 if is_MC else 0                         #True only for MC 0:disable, 1:stable only with pt>"TruthMinPT", 1+2-all, 4-single muons+tau, 8-(muon,muon)+(tau,tau) pairs
TrigRatesAlg.TruthMinPT              =300                       #cut in MeV above which to store truth particles (only stable and primary particles stored)
TrigRatesAlg.ApplyMuonCalibrations   =True                      #True for Data and MC (see link above where the tool is initialized)
TrigRatesAlg.IsRun3                  =is_Run3                   #True only for Run-3 Data 
if (do_hi2023 or do_hi2024 or do_pp2024):
  TrigRatesAlg.HLTMuonsKey           ="HLT_MuonsCB_RoI"
else:
  TrigRatesAlg.HLTMuonsKey           =""

algSeq += TrigRatesAlg
#------------------------------------------------------------


svcMgr += CfgMgr.THistSvc()
svcMgr.THistSvc.Output += ["MYSTREAM DATAFILE='myfile.root' OPT='RECREATE'"]

