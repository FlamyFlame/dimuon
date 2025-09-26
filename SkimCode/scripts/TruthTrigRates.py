InputFile = "DAOD_TRUTH0.DAOD_TRUTH0.pool.root"
m_EvtMax        =1000

import AthenaPoolCnvSvc.ReadAthenaPool
import glob
svcMgr.EventSelector.InputCollections = glob.glob(InputFile)
theApp.EvtMax = m_EvtMax


# from AthenaCommon.GlobalFlags import globalflags
# globalflags.DataSource.set_Value_and_Lock(dataSource) #data geant4
# svcMgr += CfgMgr.AthenaEventLoopMgr(EventPrintoutInterval=1)
# globalflags.DatabaseInstance = 'CONDBR2' # e.g. if you want run 2 data,  set to COMP200 if you want run1. This is completely ignored in MC.
# #globalflags.DatabaseInstance = 'COMP200' # e.g. if you want run 2 data,  set to COMP200 if you want run1. This is completely ignored in MC.
# #TODO:: WHAT ABOUT RUN3??

##------------------------------------------------------------

#the main algorithm sequence
algSeq = CfgMgr.AthSequencer("AthAlgSeq")


#------------------------------------------------------------
#main algorithm to run
TrigRatesAlg=CfgMgr.TrigRates("TrigRatesAlg",OutputLevel=INFO)
#------------------------------------------------------------



TrigRatesAlg.IsEvgen               = True
TrigRatesAlg.StoreTruth            = 1 + 2 + 4 + 8
TrigRatesAlg.StoreTruthVtx         =True
TrigRatesAlg.StoreEventInfo        = 0
TrigRatesAlg.StoreTracks           = 0
TrigRatesAlg.StorePixTracks        = 0
TrigRatesAlg.UseTrigger            = False
TrigRatesAlg.StoreVtx              = False
TrigRatesAlg.StoreL1TE             = False

algSeq += TrigRatesAlg
#------------------------------------------------------------


svcMgr += CfgMgr.THistSvc()
svcMgr.THistSvc.Output += ["MYSTREAM DATAFILE='myfile.root' OPT='RECREATE'"]

