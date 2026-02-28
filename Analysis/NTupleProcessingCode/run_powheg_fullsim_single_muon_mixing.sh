#!/bin/bash
set -eo pipefail

# Usage:
#   ./run_powheg_fullsim_single_muon_mixing.sh [bb|cc|both] [target_pairs] [run_year] [input_file] [output_file]
# Examples:
#   ./run_powheg_fullsim_single_muon_mixing.sh bb 5000 17
#   ./run_powheg_fullsim_single_muon_mixing.sh both 5000000 17

mode="${1:-both}"
target_pairs="${2:-5000000}"
run_year="${3:-17}"
input_override="${4:-}"
output_override="${5:-}"

# Prevent atlasLocalSetup from interpreting this script's positional args.
set --

cd "$PWD"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
set +e
set +o pipefail
set +u
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh"
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"
set -e
set -o pipefail
set -u

run_one_mode () {
  local mmode="$1"
  local in_file="$input_override"
  local out_file="$output_override"

  root -b -l << EOF
    .L mix_powheg_single_muon_pairs.C+
    run_powheg_single_muon_pair_mixing("${mmode}", ${target_pairs}, 1, ${run_year}, "${in_file}", "${out_file}", 0, 5.0);
    .q
EOF
}

if [[ "$mode" == "bb" || "$mode" == "cc" ]]; then
  run_one_mode "$mode"
elif [[ "$mode" == "both" ]]; then
  run_one_mode "bb"
  run_one_mode "cc"
else
  echo "Invalid mode: $mode (must be bb, cc, or both)"
  exit 1
fi
