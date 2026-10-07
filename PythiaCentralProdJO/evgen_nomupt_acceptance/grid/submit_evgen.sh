#!/bin/bash
# Grid evgen of the no-muon-pT-filter signal-acceptance sample (tracking doc
# signal_acceptance_mass_window.md D5, docs/private_pythia_sample_notes.md).
# One prun task per pT-hat slice: N accepted events/job x nJobs. Config table: grid/slices.txt
#   <cfg dir> <evts/job> <nJobs> <seedBase>
# Run from evgen_nomupt_acceptance/ :  ./grid/submit_evgen.sh [--noSubmit]
# Uses the AthGeneration 23.6.68 acm build of evgen_high_pT_additional (stock filters + an unused
# RequireSameCharge property; this JO does not set it) via --useAthenaPackages.
EXTRA="$*"
VER=v1
WD=$(pwd)
BUILD=/usatlas/u/yuhanguo/workarea/dimuon_codes/PythiaCentralProdJO/evgen_high_pT_additional
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
set --
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh > /dev/null
cd $BUILD/build && acmSetup --sourcedir=../source/athena AthGeneration,23.6.68 > /dev/null || { echo "acmSetup FAILED"; exit 1; }
lsetup panda > /dev/null
while read cfg N nJ seed; do
  [[ -z "$cfg" || "$cfg" == \#* ]] && continue
  jo=$(ls $WD/run/$cfg/mc.*.py); phys=$(basename $jo .py); phys=${phys#mc.}
  st=$WD/grid/stage/$cfg; mkdir -p $st
  cp $jo $st/; echo "evgenConfig.nEventsPerJob = $N" >> $st/$(basename $jo)
  cp -L $WD/run/Py8Control_hQCD_OSdimu_nomupt.py $st/; cp $WD/grid/run_gen.sh $st/
  outDS=user.yuhang.nomuptAcc.${phys}.${VER}
  echo "== $cfg : $N evts/job x $nJ jobs, seeds from $seed -> $outDS"
  ( cd $st && prun --exec "./run_gen.sh %RNDM=$seed $N" --nJobs $nJ --useAthenaPackages \
        --outDS $outDS --outputs EVNT.pool.root,log.generate $EXTRA < /dev/null )
done < $WD/grid/slices.txt
