#!/bin/bash
set -euo pipefail

# Condor passes $(Process) as the first argument.
process_id="${1:-0}"
file_batch=$(( process_id + 1 ))

# Optional overrides for local quick tests:
#   TARGET_PAIRS=1000 RUN_YEAR=17 ./run_powheg_fullsim_single_muon_mixing.sh 0
target_pairs="${TARGET_PAIRS:-125000}"
run_year="${RUN_YEAR:-17}"

cd "$PWD"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
set +e
set +u
set +o pipefail
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh"
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"
set -e
set -u
set -o pipefail

root -b -l << EOF
    .L mix_powheg_single_muon_pairs.C
    run_powheg_single_muon_pair_mixing("bb", ${target_pairs}, ${file_batch}, ${run_year},
        "",
        "",
        ${file_batch});
    .q
EOF
