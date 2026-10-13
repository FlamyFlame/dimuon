#!/bin/bash
# Grid evgen test for the 6 configs: one prun task per config, N events/job x nJobs.
# Run from evgen_high_pT_additional/ :  ./grid/submit_evgen.sh [--noSubmit]
# Needs a valid VOMS proxy. Uses the locally compiled GeneratorFilters patch (--useAthenaPackages).

EXTRA="$*"
VER=v1
WD=$(pwd)
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
set --   # sourced setup scripts must not see this script's arguments
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh > /dev/null
cd $WD/build && acmSetup --sourcedir=../source/athena AthGeneration,23.6.68 > /dev/null && acm compile > /dev/null || { echo "acm compile FAILED"; exit 1; }; cd $WD
lsetup panda > /dev/null
#      config dir                  evts/job nJobs seedBase  (total: 1000 / 1000 / 1000 / 1000 / 100 / 100)
CFGS=("run_sign_indep_pT125_300    250  4  120000"
      "run_same_sign_pT125_300     200  5  160000"
      "run_sign_indep_pT70_125     100 10  110000"
      "run_same_sign_pT70_125       50 20  150000"
      "run_same_sign_pT40_70         5 20  140000"
      "run_same_sign_pT24_40        10 10  130000")
for line in "${CFGS[@]}"; do
  read cfg N nJ seed <<< "$line"
  jo=$(ls $WD/run/$cfg/mc.*.py); phys=$(basename $jo .py); phys=${phys#mc.}
  st=$WD/grid/stage/$cfg; mkdir -p $st; rm -f $st/*
  cp $jo $st/; cp -L $WD/run/Py8Control_hQCD_dimu_pdf_dilepton_filters.py $st/; cp $WD/grid/run_gen.sh $st/
  outDS=user.yuhang.JOtest.${phys}.${VER}
  echo "== $cfg : $N evts/job x $nJ jobs, seeds from $seed -> $outDS"
  ( cd $st && prun --exec "./run_gen.sh %RNDM=$seed $N" --nJobs $nJ --useAthenaPackages \
        --outDS $outDS --outputs EVNT.pool.root,log.generate $EXTRA )
done
