#!/bin/bash
# Condor worker: Gen_tf <cfg> <N accepted> <seed>  -> condor_check/<cfg>/s<seed>/
cfg=$1; N=$2; seed=$3; B=/usatlas/u/yuhanguo/workarea/dimuon_codes/PythiaCentralProdJO/evgen_nomupt_acceptance
out=$B/condor_check/$cfg/s$seed; mkdir -p $out; cd $out
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
set --; source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh > /dev/null
asetup AthGeneration,23.6.68 > /dev/null
cp $B/run/$cfg/mc.*.py .; echo "evgenConfig.nEventsPerJob = $N" >> mc.*.py; cp -L $B/run/Py8Control_hQCD_OSdimu_nomupt.py .
start=$(date +%s)
Gen_tf.py --ecmEnergy=5360. --randomSeed=$seed --firstEvent=$(( seed*1000 + 1 )) --maxEvents=$N --jobConfig=$PWD --outputEVNTFile=EVNT.pool.root > gen_tf.out 2>&1
echo "EXIT $? WALL $(( $(date +%s)-start ))" >> gen_tf.out
