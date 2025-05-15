## Running test job

### Generate EVNT file
* write your python job option file; make sure file name is in the `mc.job_option_name.py` format
* set up AthGeneration release: `asetup AthGeneration,23.6.23` (or the latest)
* Run Gen_tf.py using current directory + preExec for setting commands; output EVNT file

```
Gen_tf.py --firstEvent=1 --maxEvents=10 --randomSeed=1235 \
          --jobConfig=$PWD \
          --outputEVNTFile=Pythia.EVNT.pool.root \
          --preExec "runArgs.beam1=\"PROTON\"; runArgs.beam2=\"NEUTRON\"; runArgs.ecmEnergy=5020.; runArgs.pTHatMin=60; runArgs.pTHatMax=120"
```

* Check python output messages; especially, check the configuration is as expected + that the run arguments are correctly set
* use `checkFile.py` to check EVNT file
  * PDF weights

### Check output: ENVT --> AOD --> NTUP
* EVNT --> AOD: use Reco_tf.py

```
Derivation_tf.py \
  --inputEVNTFile Pythia.EVNT.pool.root \
  --outputDAODFile DAOD_TRUTH0.pool.root \
  --formats TRUTH0 \
  --maxEvents 1000
```

* use `checkxAOD.py` to check content of the AOD file
* From the AOD, access the event weight using the `EventInfo::mcEventWeights()` function. (This is already done in the `TrigRates` skimming-code package, and the info is saved as a branch `EventWeights`.)
See [AnalysisWorkBook](https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/PhysicsAnalysisWorkBookRel20MC#Accessing_generator_event_weight)

* Run skimming code to get ntuple from AOD
  * Set up the analysis package

```
cd /afs/cern.ch/user/y/yuhang/eos/dimuon/SkimCode
source setup_25.sh
```

  * Run python configuration file: `athena TrigRates.py`

* Compare output ntuple with ntuple from private sample with same collision energy & pTHat range
