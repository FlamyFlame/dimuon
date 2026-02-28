#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"

# Setup ATLAS environment
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh"

# Setup LCG
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

input_file="/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp17.bb.Feb2026.v1._MYSTREAM/single_muon_trees_powheg_bb_fullsim_pp17.root"
output_file="/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample/user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp17.bb.Feb2026.v1._MYSTREAM/mixed/muon_pairs_powheg_bb_fullsim_mixed_batch999_test.root"

# Tiny test target: total 30 mixed pairs (20 with truth_pair_pt<=20, 10 with truth_pair_pt>20)
root -b -l << EOF
    .L mix_powheg_single_muon_pairs.C
    run_powheg_single_muon_pair_mixing("bb", 20, 10, 999, 17,
      "${input_file}",
      "${output_file}",
      999);
    .q
EOF

echo "Tiny bb test output: ${output_file}"
