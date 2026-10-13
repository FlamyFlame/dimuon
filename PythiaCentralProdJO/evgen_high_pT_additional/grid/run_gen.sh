#!/bin/bash
# Grid worker: Gen_tf for N accepted events with the single mc.*.py in this directory.
# Args: <seed> <nEvents>
seed=$1; N=$2
start=$(date +%s)
Gen_tf.py --ecmEnergy=5360. --randomSeed=$seed --firstEvent=$(( seed*1000 + 1 )) --maxEvents=$N \
          --jobConfig=$PWD --outputEVNTFile=EVNT.pool.root
rc=$?
echo "Job completed in $(( $(date +%s)-start )) seconds, Gen_tf rc=$rc" | tee -a log.generate
exit $rc
