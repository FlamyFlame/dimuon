#!/bin/bash

# The argument gets passed as the argument $(Process) in batch-submission script
file_batch=$(( $1 + 1 ))

cd $PWD

# Setup ATLAS environment
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

# Setup LCG
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

# Run bb mixing with fixed target: 50k mixed pairs per job (truth-pt^5 weighted sampling)
root -b -l << EOF
    .L mix_powheg_single_muon_pairs.C
  run_powheg_single_muon_pair_mixing("bb", 50000, $file_batch, 17,
      "/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp17.bb.Feb2026.v1._MYSTREAM/single_muon_trees_powheg_bb_fullsim_pp17.root",
      "",
      $file_batch,
      5.0);
    .q
EOF
