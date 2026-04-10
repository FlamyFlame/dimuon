#!/usr/bin/env python

import argparse
import os

from AthenaConfiguration.AllConfigFlags import initConfigFlags
from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.MainServicesConfig import MainServicesCfg
from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg

from MuonSelectorTools.MuonSelectorToolsConfig import MuonSelectionToolCfg
from MuonMomentumCorrections.MCastCfg import setupMCastToolCfg
from MuonEfficiencyCorrections.MuonEfficiencyCorrectionsCfg import MuonEfficiencyCorrectionsCfg
from TrigDecisionTool.TrigDecisionToolConfig import TrigDecisionToolCfg


# data with muons
do_pp2015 = False
do_pp2017 = False
do_pp2024 = False
do_hi2018 = False
do_hi2015 = False
do_hi2023 = False
do_hi2024 = False
do_hi2025 = False

# overlay
do_pp_MC_fullsim_17 = False
do_pp_MC_fullsim_24 = False

# Resolve run mode (priority order):
#   1) explicit TRIGRATES_RUNMODE environment variable (most robust for grid)
#   2) infer from --filesInput dataset name
#   3) infer from run directory name (local convenience)
_cli_parser = argparse.ArgumentParser(add_help=False)
_cli_parser.add_argument("--filesInput", type=str, default="")
_cli_args, _ = _cli_parser.parse_known_args()
_files_input_arg = (_cli_args.filesInput or "").lower()
_run_mode_arg = (os.environ.get("TRIGRATES_RUNMODE", "") or "").strip().lower()


def _set_run_mode(mode):
	global do_hi2023, do_hi2024, do_hi2025, do_pp2024, do_pp_MC_fullsim_24
	do_hi2023 = (mode == "hi2023")
	do_hi2024 = (mode == "hi2024")
	do_hi2025 = (mode == "hi2025")
	do_pp2024 = (mode == "pp2024")
	do_pp_MC_fullsim_24 = (mode == "ppmcfullsim2024")


if _run_mode_arg in ("hi2023", "hi2024", "hi2025", "pp2024", "ppmcfullsim2024"):
	_set_run_mode(_run_mode_arg)
elif _files_input_arg:
	if "data23_hi" in _files_input_arg:
		_set_run_mode("hi2023")
	elif "data24_hi" in _files_input_arg:
		_set_run_mode("hi2024")
	elif "data25_hi" in _files_input_arg:
		_set_run_mode("hi2025")
	elif "data24_5p36tev" in _files_input_arg or "data24_pp" in _files_input_arg:
		_set_run_mode("pp2024")
else:
	_run_dir = os.path.basename(os.getcwd())
	if _run_dir == "run_23hi":
		_set_run_mode("hi2023")
	elif _run_dir == "run_24hi":
		_set_run_mode("hi2024")
	elif _run_dir == "run_25hi":
		_set_run_mode("hi2025")
	elif _run_dir == "run_24pp":
		_set_run_mode("pp2024")


# ------------------------------------------------------------
# whether the sample is HION or not
is_HION = False
if (do_hi2018 or
	 do_hi2015 or
	 do_hi2023 or
	 do_hi2024 or
	 do_hi2025
	 ):
	is_HION = True
# ------------------------------------------------------------


# ------------------------------------------------------------
# for Run3 data
is_Run3 = False
if (do_hi2023 or do_hi2024 or do_hi2025 or do_pp2024 or do_pp_MC_fullsim_24):
	is_Run3 = True
# ------------------------------------------------------------

# ------------------------------------------------------------
is_MC = do_pp_MC_fullsim_17 or do_pp_MC_fullsim_24
# ------------------------------------------------------------


# CA config here is intended for modern Run-3-style configurations
if not (do_hi2023 or do_hi2024 or do_hi2025 or do_pp2024 or do_pp_MC_fullsim_24):
	print("*"*50, "\nTrigRates_CA.py is for modern Run-3 style samples. Use TrigRates_JO.py for older samples.\n", "*"*50)
	exit()


RunYear = 0
m_EvtMax = 1000
dataSource = 'data'
GRL = []
MinBias_triggers = []
Muon_triggers = []
DiMuon_triggers = []

if do_hi2023:
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
elif do_hi2024:
	RunYear = 2024
	GRL = ["physics_HI2024_50ns.xml"]
	InputFile = "/eos/user/y/yuhang/data/data24_hi_testfile_AOD/data24_hi.00488980.physics_HardProbes.merge.AOD.f1537_m2259._lb0230._0001.1"
	Muon_triggers = ["HLT_mu4_L1MU3V",
									 "HLT_mu6_L1MU3V",
									 "HLT_mu6_L1MU5VF",
									 "HLT_mu8_L1MU5VF",
									 "HLT_mu10_L1MU8F",
									 "HLT_mu10_L1MU5VF",
									 "HLT_mu4noL1_hi_uccTh3_L1jTE10000",
									 "HLT_mu4noL1_hi_uccTh2_L1jTE9000",
									 "HLT_mu4noL1_L1ZDC_HELT25_jTE4000",
									 "HLT_mu4noL1_L1ZDC_HELT20_jTE4000",
									 "HLT_mu4noL1_L1ZDC_HELT15_jTE4000"]
	DiMuon_triggers = ["HLT_2mu4_L12MU3V", "HLT_mu4_mu4noL1_L1MU3V"]
elif do_hi2025:
	RunYear = 2025
	GRL = ["physics_HI2025_50ns_PbPb_IgnoreBSPOT_INVALID.xml"]
	# data25_hi Pb+Pb HardProbes AOD (periods R/S/T, runs 510493+).
	# Files are not yet staged at CERN EOS; update path once staged or use xrootd/Rucio.
	InputFile = "/eos/atlas/atlastier0/rucio/data25_hi/physics_HardProbes/00510493/data25_hi.00510493.physics_HardProbes.merge.AOD.f1655_m2272/data25_hi.00510493.physics_HardProbes.merge.AOD.f1655_m2272._lb0001._0001.1"
	Muon_triggers = ["HLT_mu4_L1MU3V",
									 "HLT_mu6_L1MU3V",
									 "HLT_mu6_L1MU5VF",
									 "HLT_mu8_L1MU5VF",
									 "HLT_mu10_L1MU8F",
									 "HLT_mu10_L1MU5VF"]
	DiMuon_triggers = ["HLT_2mu4_L12MU3V", "HLT_mu4_mu4noL1_L1MU3V"]
elif do_pp2024:
	GRL = ["physics_2024ppRef_25ns.xml"]
	InputFile = "/eos/user/y/yuhang/data/data24_pp_testfile_AOD/data24_5p36TeV.00488427.physics_Main.merge.AOD.f1529_m2259._lb0100._0004.1"
	Muon_triggers = ["HLT_mu4_L1MU3V",
									 "HLT_mu6_L1MU3V",
									 "HLT_mu8_L1MU5VF",
									 "HLT_mu10_L1MU8F",
									 "HLT_mu12_L1MU8F",
									 "HLT_mu15_L1MU8F",
									 "HLT_mu15_L1MU14FCH"]
	DiMuon_triggers = ["HLT_2mu3_L12MU3V", "HLT_2mu4_L12MU3V", "HLT_mu4_mu6_L12MU3V", "HLT_mu4_mu4noL1_L1MU3V"]
elif do_pp_MC_fullsim_24:
	dataSource = 'geant4'
	InputFile = "/eos/user/y/yuhang/data/mc_powheg_fullsim/example_fullsim_AOD_files/AOD.36949837._016835.pool.root.1"
else:
	print("*"*50, "\nUnknown Dataset\n", "*"*50)
	exit()

print(MinBias_triggers + Muon_triggers + DiMuon_triggers)


def build_cfg(evt_max=1000):
	flags = initConfigFlags()
	flags.Input.Files = [InputFile]   # local default; overridden by --filesInput=%IN on the grid
	flags.Input.isMC = is_MC
	flags.Exec.MaxEvents = evt_max    # local default; overridden by --evtMax=%MAXEVENTS on the grid
	# Pick up --filesInput / --evtMax / other CA flags passed on the command line.
	# This is required when submitted via:
	#   pathena --trf "TRIGRATES_RUNMODE=hi2023 athena.py TrigRates_CA.py --filesInput=%IN --evtMax=%MAXEVENTS" ...
	flags.fillFromArgs()
	flags.lock()

	cfg = MainServicesCfg(flags)
	cfg.merge(PoolReadCfg(flags))

	cfg.addService(CompFactory.THistSvc(Output=["MYSTREAM DATAFILE='myfile.root' OPT='RECREATE'"]))

	trk_sel_cls = CompFactory.getComp("InDet::InDetTrackSelectionTool")
	grl_cls = CompFactory.getComp("GoodRunsListSelectionTool")
	match_cls = CompFactory.getComp("Trig::R3MatchingTool")

	trk_minbias = trk_sel_cls("TrackSelectionTool_MinBias", CutLevel="MinBias", minPt=100.)
	trk_hiloose = trk_sel_cls("TrackSelectionTool_HILoose", CutLevel="HILoose", minPt=100.)
	trk_hitight = trk_sel_cls("TrackSelectionTool_HITight", CutLevel="HITight", minPt=100.)
	cfg.addPublicTool(trk_minbias)
	cfg.addPublicTool(trk_hiloose)
	cfg.addPublicTool(trk_hitight)

	mu_sel = cfg.popToolsAndMerge(MuonSelectionToolCfg(flags, name="MyMuonSelectionTool", MaxEta=2.5, IsRun3Geo=is_Run3))
	if is_HION:
		mu_sel.TrtCutOff = True
	cfg.addPublicTool(mu_sel)

	mu_calib = cfg.popToolsAndMerge(setupMCastToolCfg(flags, name="MyMuonCalibrationAndSmearingTool", calibMode=0, IsRun3Geo=is_Run3))
	cfg.addPublicTool(mu_calib)

	eff_release_run3 = "250418_Preliminary_r24run3"
	eff_medium = cfg.popToolsAndMerge(MuonEfficiencyCorrectionsCfg(
		flags,
		name="MuonEfficiencyTool_Medium",
		WorkingPoint="Medium",
		CalibrationRelease=eff_release_run3,
	))
	eff_tight = cfg.popToolsAndMerge(MuonEfficiencyCorrectionsCfg(
		flags,
		name="MuonEfficiencyTool_Tight",
		WorkingPoint="Tight",
		CalibrationRelease=eff_release_run3,
	))
	cfg.addPublicTool(eff_medium)
	cfg.addPublicTool(eff_tight)

	grl_tool = grl_cls("MyGRLTool", GoodRunsListVec=GRL)
	cfg.addPublicTool(grl_tool)

	trig_match_tool = None
	if not is_MC:
		tdt = cfg.getPrimaryAndMerge(TrigDecisionToolCfg(flags))
		trig_match_tool = match_cls("MyTriggerMatchTool", TrigDecisionTool=tdt)
		cfg.addPublicTool(trig_match_tool)

	import HFtrigValidation.HFtrigValidationConf  # noqa: F401
	trigrates_cls = CompFactory.getComp("TrigRates")
	alg = trigrates_cls("TrigRatesAlg", OutputLevel=3)

	alg.TrackSelectionTool_MinBias = trk_minbias
	alg.TrackSelectionTool_HILoose = trk_hiloose
	alg.TrackSelectionTool_HITight = trk_hitight
	alg.MuonSelectionTool = mu_sel
	if not is_MC:
		alg.TriggerMatchTool = trig_match_tool

	alg.GRLTool = grl_tool
	alg.muonCorrectionTool = mu_calib
	alg.use_effi_SF_tool = True
	alg.effi_SF_tool_medium = eff_medium
	alg.effi_SF_tool_tight = eff_tight

	alg.RunYear = RunYear
	alg.IsEvgen = False
	alg.UseGRL = True
	alg.StoreAllEvents = False
	alg.UseTrigger = (not is_MC)
	alg.StoreL1Decision = (not is_MC)
	alg.TriggerChains = "|".join(MinBias_triggers)
	alg.MuonTriggerChains = "|".join(Muon_triggers)
	alg.DiMuonTriggerChains = "|".join(DiMuon_triggers)
	alg.StoreEventInfo = 1
	# StoreTracks bitmap: 1=trk_numqual only (all+PPMinBias+HILoose+HITight counts, with and without pt>400 cut)
	#                     1+2=also per-track vectors (pt/eta/phi/charge/quality)
	#                     1+2+4=also hit details (d0, z0, pixel/SCT hits, etc.)
	#                     1+2+4+8=also truth links (MC only)
	alg.StoreTracks = 1
	alg.StorePixTracks = 0
	alg.MaxZvtx = 250
	alg.StoreVtx = True
	alg.StoreL1TE = (not is_MC)
	alg.StoreMuonTruth = is_MC
	alg.StoreSingleMuon = True
	alg.StoreAcoplanarMuon = True
	alg.HIEventShapeContainerKey = "HIEventShape"
	alg.METContainerKey = ""
	alg.StoreTruthVtx = is_MC
	alg.StoreTruth = 1 + 2 + 4 + 8 if is_MC else 0
	alg.TruthMinPT = 300
	alg.ApplyMuonCalibrations = True
	alg.IsRun3 = is_Run3
	alg.IsZdcCalib = False
	alg.StoreZdc = 1 if is_HION else 0   # 1=basic ZDC (energy/time/status/PreSampleAmp); add 2 for RPD centroid data
	alg.ZdcAuxSuffix = ""
	if do_hi2023 or do_hi2024 or do_hi2025 or do_pp2024:
		alg.HLTMuonsKey = "HLT_MuonsCB_RoI"
	else:
		alg.HLTMuonsKey = ""

	cfg.addEventAlgo(alg)
	return cfg


def _parse_evtmax(default_evtmax):
	parser = argparse.ArgumentParser(add_help=False)
	parser.add_argument("--evtMax", type=int, default=default_evtmax)
	parser.add_argument("-n", "--maxEvents", type=int, default=None)
	args, _ = parser.parse_known_args()
	if args.maxEvents is not None:
		return args.maxEvents
	return args.evtMax


if __name__ == "__main__":
	import sys
	from AthenaCommon.Configurable import ConfigurableCABehavior

	evtmax = _parse_evtmax(m_EvtMax)

	# On the grid the PanDA pilot pre-loads legacy Athena modules
	# (AthenaCommon/Preparation.py, Execution.py, FakeAppMgr.py) before our
	# script runs.  These cause two problems inside ConfigurableCABehavior:
	#
	# Problem 1 – stale AthSequencer type in CFElements:
	#   CFElements.py executes  AthSequencer = CompFactory.AthSequencer  at
	#   import time, while still in legacy mode.  CompFactory.__getattr__ is
	#   mode-aware but the cached name is not, so isSequence() keeps checking
	#   isinstance(obj, <legacy class>).  When MainServicesCfg then creates
	#   CompFactory.AthSequencer('AthAlgEvtSeq') in CA mode it produces a
	#   GaudiConfig2 instance that fails the legacy isinstance check →
	#   "AthAlgEvtSeq is not a sequence".
	#   Fix: re-bind CFElements.AthSequencer to the CA class while inside the
	#   CA context (where CompFactory.__getattr__ already returns the CA class).
	#
	# Problem 2 – legacy instances in allConfigurables:
	#   FakeAppMgr.py creates legacy AthAlgEvtSeq / AthMasterSeq instances.
	#   Fix: clear the registry so MainServicesCfg builds them fresh in CA mode.
	with ConfigurableCABehavior():
		# Fix problem 1: update the stale CFElements cache
		import AthenaCommon.CFElements as _cfe
		from AthenaConfiguration.ComponentFactory import CompFactory as _cf
		_cfe.AthSequencer = _cf.AthSequencer  # now resolves to GaudiConfig2 class

		# Fix problem 2: discard legacy instances from the grid preamble
		from AthenaCommon.Configurable import Configurable
		Configurable.allConfigurables.clear()

		cfg = build_cfg(evtmax)
		sc = cfg.run()

	sys.exit(not sc.isSuccess())
