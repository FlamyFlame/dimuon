#!/bin/bash
# local Gen_tf test: local_test/run_local.sh <cfg> <nEvents> <seed>
cfg=$1; N=$2; seed=$3; B=/usatlas/u/yuhanguo/workarea/dimuon_codes/PythiaCentralProdJO/evgen_nomupt_acceptance
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
set --; source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh > /dev/null
out=$B/local_test/${cfg}_s$seed; mkdir -p $out; cd $out
asetup AthGeneration,23.6.68 > /dev/null
out=$B/local_test/${cfg}_s$seed; mkdir -p $out; cd $out
start=$(date +%s)
Gen_tf.py --ecmEnergy=5360. --randomSeed=$seed --firstEvent=1 --maxEvents=$N --jobConfig=$B/run/$cfg --outputEVNTFile=test.EVNT.pool.root > gen_tf.out 2>&1
echo "EXIT $? WALL $(( $(date +%s)-start ))" >> gen_tf.out
