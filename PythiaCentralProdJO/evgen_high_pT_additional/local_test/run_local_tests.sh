#!/bin/bash
# Local Gen_tf test: N accepted events per configuration, all 6 configs in parallel.
# Usage: local_test/run_local_tests.sh [nEvents=10] [config dir names...]  (default: all run/run_*)
N=${1:-10}; shift
CFGS="$*"
JOTEST_WD=$(cd "$(dirname "$0")/.." && pwd)
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh > /dev/null
cd $JOTEST_WD/build && acmSetup --sourcedir=../source/athena AthGeneration,23.6.68 > /dev/null && cd $JOTEST_WD
seed=1000
for jo in $( [ -n "$CFGS" ] && for c in $CFGS; do echo $JOTEST_WD/run/$c; done || ls -d $JOTEST_WD/run/run_* ); do
  cfg=$(basename $jo); seed=$((seed+1)); [ -n "$CFGS" ] && seed=$((seed+100))
  out=$JOTEST_WD/local_test/$cfg; mkdir -p $out
  ( cd $out && start=$(date +%s) && \
    Gen_tf.py --ecmEnergy=5360. --randomSeed=$seed --firstEvent=1 --maxEvents=$N \
              --jobConfig=$jo --outputEVNTFile=test.EVNT.pool.root > gen_tf.out 2>&1
    rc=$?; echo "EXIT $rc WALL $(( $(date +%s)-start ))" >> gen_tf.out ) &
done
wait
echo ALL_DONE
