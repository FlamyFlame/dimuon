# Legacy JO configuration for TrigRates (pre-R24)
import os

# data with muons
do_pp2015 = False
do_pp2017 = False
do_pp2024 = False
do_hi2018 = False
do_hi2015 = False
do_hi2023 = False
do_hi2024 = True

# overlay
do_pp_MC_fullsim_17 = False
do_pp_MC_fullsim_24 = False

_run_dir = os.path.basename(os.getcwd())
if _run_dir == "run_powheg_pp17fullsim":
	do_pp2015 = False
	do_pp2017 = False
	do_pp2024 = False
	do_hi2018 = False
	do_hi2015 = False
	do_hi2023 = False
	do_hi2024 = False
	do_pp_MC_fullsim_17 = True
	do_pp_MC_fullsim_24 = False


# ------------------------------------------------------------
# whether the sample is HION or not
is_HION = False
if (do_hi2018 or
	 do_hi2015 or
	 do_hi2023 or
	 do_hi2024
	 ):
	is_HION = True
# ------------------------------------------------------------


# ------------------------------------------------------------
# for Run3 data
is_Run3 = False
if (do_hi2023 or do_hi2024 or do_pp2024 or do_pp_MC_fullsim_24):
	is_Run3 = True
# ------------------------------------------------------------


# ------------------------------------------------------------
# for Analysis in Releases>=22 and (athena) Releases>=24
is_Rel22orAbove = False
if (do_hi2023 or do_hi2024 or do_pp2024 or do_pp_MC_fullsim_24):
	is_Rel22orAbove = True
# ------------------------------------------------------------

# ------------------------------------------------------------
is_MC = do_pp_MC_fullsim_17 or do_pp_MC_fullsim_24
# ------------------------------------------------------------

# Legacy JO is intended for pre-R24 only
if do_hi2024 or do_pp2024 or do_pp_MC_fullsim_24:
	print("*"*50, "\nTrigRates_JO.py is for pre-R24 only. Use TrigRates_CA.py for R24/R25.\n", "*"*50)
	exit()


# ------------------------------------------------------------
RunYear = 0
m_EvtMax = 1000
dataSource = 'data'
GRL = []
MinBias_triggers = []
Muon_triggers = []
DiMuon_triggers = []

if do_hi2018:
	RunYear = 2018
	GRL = ["data18_hi.periodAllYear_DetStatus-v106-pro22-14_Unknown_PHYS_HeavyIonP_All_Good.xml"]
	InputFile = "/afs/cern.ch/user/s/soumya/workarea/DATA/JobTestData/data18_hi/data18_hi.00366142.physics_HardProbes.merge.AOD.f1027_m2037._lb0570._0003.1"
	MinBias_triggers = []
	DiMuon_triggers = ["HLT_2mu3", "HLT_2mu4", "HLT_2mu6", "HLT_2mu8", "HLT_mu4_mu4noL1"]
elif do_hi2015:
	RunYear = 2015
	GRL = ["data15_hi.periodAllYear_DetStatus-v75-repro20-01_DQDefects-00-02-02_PHYS_HeavyIonP_All_Good.xml"]
	InputFile = "/afs/cern.ch/user/s/soumya/workarea/DATA/JobTestData/data15_hi/AOD.16615737._000065.pool.root.1"
	MinBias_triggers = []
	DiMuon_triggers = ["HLT_2mu3", "HLT_2mu4", "HLT_2mu6", "HLT_2mu8", "HLT_mu4_mu4noL1"]
elif do_hi2023:
	RunYear = 2023
	GRL = ["data23_hi.periodAllYear_DetStatus-v113-pro31-08_MERGED_PHYS_HeavyIonP_All_Good.xml"]
	InputFile = "/eos/user/y/yuhang/data/data23_hi_testfile_AOD/data23_hi.00462240.physics_HardProbes.merge.AOD.f1399_m2209._lb0409._0002.1"
	Muon_triggers = ["HLT_mu4_L1MU3V",
									 "HLT_mu6_L1MU3V",
									 "HLT_mu6_L1MU5VF",
									 "HLT_mu4_L1MU3V_VTE50",
									 "HLT_mu6_L1MU3V_VTE50",
									 "HLT_mu8_L1MU5VF_VTE50"]
	DiMuon_triggers = ["HLT_2mu4_L12MU3V", "HLT_mu4_mu4noL1_L1MU3V"]
elif do_pp2017:
	GRL = ["data17_5TeV.periodAllYear_DetStatus-v98-pro21-16_Unknown_PHYS_StandardGRL_All_Good_25ns_ignore_GLOBAL_LOWMU.xml"]
	InputFile = "/eos/user/y/yuhang/data/pp_17/data17_5TeV/*"
	Muon_triggers = ["HLT_mu4", "HLT_mu6", "HLT_mu8", "HLT_mu10"]
	DiMuon_triggers = ["HLT_2mu4"]
elif do_pp2015:
	GRL = ["data15_5TeV.periodAllYear_DetStatus-v75-repro20-01_DQDefects-00-02-02_PHYS_HeavyIonP_All_Good.xml",
				 "data15_5TeV.periodVdM_DetStatus-v75-repro20-01_DQDefects-00-02-02_PHYS_HeavyIonP_All_Good.xml"]
	InputFile = "/eos/user/y/yuhang/data/pp_15/data15_5TeV_r9582/*"
	Muon_triggers = ["HLT_mu4", "HLT_mu6", "HLT_mu8", "HLT_mu10"]
	DiMuon_triggers = ["HLT_2mu4"]
elif do_pp_MC_fullsim_17:
	dataSource = 'geant4'
	InputFile = "/eos/user/y/yuhang/data/mc_powheg_fullsim/example_fullsim_AOD_files/AOD.36949837._016835.pool.root.1"
else:
	print("*"*50, "\nUnknown Dataset\n", "*"*50)
	exit()
print(MinBias_triggers + Muon_triggers + DiMuon_triggers)
# ------------------------------------------------------------


import AthenaPoolCnvSvc.ReadAthenaPool
import glob
svcMgr.EventSelector.InputCollections = glob.glob(InputFile)
theApp.EvtMax = m_EvtMax


from AthenaCommon.GlobalFlags import globalflags
globalflags.DataSource.set_Value_and_Lock(dataSource)
svcMgr += CfgMgr.AthenaEventLoopMgr(EventPrintoutInterval=1)
globalflags.DatabaseInstance = 'CONDBR2'


# the tool service
from AthenaCommon.AppMgr import ToolSvc

# ------------------------------------------------------------
# Configure Track selection tool
ToolSvc += CfgMgr.InDet__InDetTrackSelectionTool("TrackSelectionTool_MinBias", CutLevel="MinBias", minPt=100.)
ToolSvc += CfgMgr.InDet__InDetTrackSelectionTool("TrackSelectionTool_HILoose", CutLevel="HILoose", minPt=100.)
ToolSvc += CfgMgr.InDet__InDetTrackSelectionTool("TrackSelectionTool_HITight", CutLevel="HITight", minPt=100.)
# ------------------------------------------------------------

# ---------------------------- Muon Selection Tool --------------------------------
if not hasattr(ToolSvc, "MyMuonSelectionTool"):
	from MuonSelectorTools.MuonSelectorToolsConf import CP__MuonSelectionTool
	ToolSvc += CP__MuonSelectionTool("MyMuonSelectionTool")
	ToolSvc.MyMuonSelectionTool.MaxEta = 2.5

	if is_HION:
		ToolSvc.MyMuonSelectionTool.TrtCutOff = True

	if is_Rel22orAbove:
		ToolSvc.MyMuonSelectionTool.IsRun3Geo = False
		if is_Run3:
			ToolSvc.MyMuonSelectionTool.IsRun3Geo = True
	print(ToolSvc.MyMuonSelectionTool)
# ------------------------------------------------------------

# ---------------------------- Muon Calibration Tool --------------------------------
if is_Rel22orAbove:
	ToolSvc += CfgMgr.CP__MuonCalibTool("MyMuonCalibrationAndSmearingTool")
	ToolSvc.MyMuonCalibrationAndSmearingTool.calibMode = 0
	ToolSvc.MyMuonCalibrationAndSmearingTool.IsRun3Geo = False
	if is_Run3:
		ToolSvc.MyMuonCalibrationAndSmearingTool.IsRun3Geo = True
else:
	from MuonMomentumCorrections.MuonMomentumCorrectionsConf import CP__MuonCalibrationAndSmearingTool
	ToolSvc += CP__MuonCalibrationAndSmearingTool("MyMuonCalibrationAndSmearingTool")
	ToolSvc.MyMuonCalibrationAndSmearingTool.Year = "Data18"
	ToolSvc.MyMuonCalibrationAndSmearingTool.StatComb = False
	ToolSvc.MyMuonCalibrationAndSmearingTool.SagittaRelease = "sagittaBiasDataAll_03_02_19_Data18"
	ToolSvc.MyMuonCalibrationAndSmearingTool.sagittaMapsInputType = 0
	ToolSvc.MyMuonCalibrationAndSmearingTool.SagittaCorr = True
	ToolSvc.MyMuonCalibrationAndSmearingTool.doSagittaMCDistortion = False
	ToolSvc.MyMuonCalibrationAndSmearingTool.Release = "Recs2020_03_03"
	ToolSvc.MyMuonCalibrationAndSmearingTool.SagittaCorrPhaseSpace = True

	if do_hi2015:
		ToolSvc.MyMuonCalibrationAndSmearingTool.Year = "Data16"
	elif do_hi2018:
		pass
	elif do_pp2017:
		ToolSvc.MyMuonCalibrationAndSmearingTool.Year = "Data17"
		ToolSvc.MyMuonCalibrationAndSmearingTool.SagittaRelease = "sagittaBiasDataAll_03_02_19_Data17"
	elif do_pp2015:
		ToolSvc.MyMuonCalibrationAndSmearingTool.Year = "Data16"
		ToolSvc.MyMuonCalibrationAndSmearingTool.SagittaRelease = "sagittaBiasDataAll_03_02_19_Data17"
	elif do_pp_MC_fullsim_17:
		ToolSvc.MyMuonCalibrationAndSmearingTool.Year = "Data17"
		ToolSvc.MyMuonCalibrationAndSmearingTool.SagittaRelease = "sagittaBiasDataAll_03_02_19_Data17"
	else:
		print("*"*50, "\nCould not configure MuonCalibrationAndSmearingTool\n", "*"*50)
		exit()

print(ToolSvc.MyMuonCalibrationAndSmearingTool)
# ------------------------------------------------------------

# ---------------------------- Muon Efficiency Tool --------------------------------
from MuonEfficiencyCorrections.CommonToolSetup import *
scale_medium = GetMuonEfficiencyTool("Medium", Release="210222_Precision_r21")
scale_tight = GetMuonEfficiencyTool("Tight", Release="210222_Precision_r21")
# ------------------------------------------------------------

# ---------------------------- Muon Trigger Matching Tool --------------------------------
if not is_MC:
	if is_Run3:
		ToolSvc += CfgMgr.Trig__R3MatchingTool("MyTriggerMatchTool", OutputLevel=DEBUG)
	else:
		ToolSvc += CfgMgr.Trig__MatchingTool("MyTriggerMatchTool", OutputLevel=DEBUG)
# ------------------------------------------------------------

# ---------------------------- GRL Tool --------------------------------
MyGRLTool = CfgMgr.GoodRunsListSelectionTool("MyGRLTool")
MyGRLTool.GoodRunsListVec = GRL
ToolSvc += MyGRLTool
# ------------------------------------------------------------

MessageSvc.defaultLimit = 1000

# the main algorithm sequence
algSeq = CfgMgr.AthSequencer("AthAlgSeq")

# ------------------------------------------------------------
# main algorithm to run
TrigRatesAlg = CfgMgr.TrigRates("TrigRatesAlg", OutputLevel=INFO)
# ------------------------------------------------------------
TrigRatesAlg.TrackSelectionTool_MinBias = ToolSvc.TrackSelectionTool_MinBias
TrigRatesAlg.TrackSelectionTool_HILoose = ToolSvc.TrackSelectionTool_HILoose
TrigRatesAlg.TrackSelectionTool_HITight = ToolSvc.TrackSelectionTool_HITight
TrigRatesAlg.MuonSelectionTool = ToolSvc.MyMuonSelectionTool
if not is_MC:
	TrigRatesAlg.TriggerMatchTool = ToolSvc.MyTriggerMatchTool

TrigRatesAlg.GRLTool = ToolSvc.MyGRLTool
TrigRatesAlg.muonCorrectionTool = ToolSvc.MyMuonCalibrationAndSmearingTool
TrigRatesAlg.use_effi_SF_tool = True
TrigRatesAlg.effi_SF_tool_medium = scale_medium
TrigRatesAlg.effi_SF_tool_tight = scale_tight

TrigRatesAlg.RunYear = RunYear
TrigRatesAlg.IsEvgen = False
TrigRatesAlg.UseGRL = True
TrigRatesAlg.StoreAllEvents = False
TrigRatesAlg.UseTrigger = (not is_MC)
TrigRatesAlg.StoreL1Decision = (not is_MC)
TrigRatesAlg.TriggerChains = "|".join(MinBias_triggers)
TrigRatesAlg.MuonTriggerChains = "|".join(Muon_triggers)
TrigRatesAlg.DiMuonTriggerChains = "|".join(DiMuon_triggers)
TrigRatesAlg.StoreEventInfo = 1
TrigRatesAlg.StoreTracks = 0
TrigRatesAlg.StorePixTracks = 0
TrigRatesAlg.MaxZvtx = 250
TrigRatesAlg.StoreVtx = True
TrigRatesAlg.StoreL1TE = (not is_MC)
TrigRatesAlg.StoreMuonTruth = is_MC
TrigRatesAlg.StoreSingleMuon = True
TrigRatesAlg.StoreAcoplanarMuon = True
TrigRatesAlg.HIEventShapeContainerKey = "HIEventShape"
TrigRatesAlg.METContainerKey = ""
TrigRatesAlg.StoreTruthVtx = is_MC
TrigRatesAlg.StoreTruth = 1 + 2 + 4 + 8 if is_MC else 0
TrigRatesAlg.TruthMinPT = 300
TrigRatesAlg.ApplyMuonCalibrations = True
TrigRatesAlg.IsRun3 = is_Run3
if do_hi2023:
	TrigRatesAlg.HLTMuonsKey = "HLT_MuonsCB_RoI"
else:
	TrigRatesAlg.HLTMuonsKey = ""

algSeq += TrigRatesAlg
# ------------------------------------------------------------

svcMgr += CfgMgr.THistSvc()
svcMgr.THistSvc.Output += ["MYSTREAM DATAFILE='myfile.root' OPT='RECREATE'"]
