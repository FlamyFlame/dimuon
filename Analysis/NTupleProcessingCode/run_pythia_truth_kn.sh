#!/bin/bash
# Run Pythia truth analysis (CRTP). Queue index maps to kinematic range / batch.
# Usage: ./run_pythia_truth_kn.sh <kn_index> [is_private] [e_com]
#   kn_index:   0-5 (kinematic range for non-private); for private batch_num = kn_index + 1 (valid: 0,1 -> batch 1,2)
#   is_private: 0 = non-private (central), 1 = private (default 0)
#   e_com:      center-of-mass energy in TeV for non-private mode (5.36 or 5.02; default 5.36)
#
# Examples:
#   ./run_pythia_truth_kn.sh 0        -> kn=0, non-private, e_com=5.36
#   ./run_pythia_truth_kn.sh 2 1      -> batch_num=3, private
#   ./run_pythia_truth_kn.sh 4 0 5.02 -> kn=4, non-private, e_com=5.02

kn=${1:-0}
is_private=${2:-0}
e_com=${3:-5.36}

cd "$(dirname "$0")"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

# batch_num: for private use kn+1 (kn=0->batch1, kn=1->batch2);
# for non-private use kn directly (constructor now sets kn_batch from batch_num)
if [[ "$is_private" -eq 1 ]]; then
	batch_num=$(( kn + 1 ))
else
	batch_num=$kn
fi

echo "Starting PythiaTruthAnalysis: kn=${kn}, is_private=${is_private}, batch_num=${batch_num}, e_com=${e_com}, nevents_max=${nevents_max}"

root -b -l << EOF
	.L PythiaAnalysisClasses.h
	PythiaTruthAnalysis pw($batch_num, (bool)$is_private, $e_com);
	pw.Run();
	.q;
EOF
