#!/bin/bash
# Run Pythia truth analysis (non-private, pp isospin only).
# Identical to run_pythia_truth_kn.sh but passes pp_only=true.
# Usage: ./run_pythia_truth_kn_pp_only.sh <kn_index> [e_com]
#   kn_index: 0-5
#   e_com:    5.36 (default) or 5.02

kn=${1:-0}
e_com=${2:-5.36}

cd "$(dirname "$0")"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

echo "Starting PythiaTruthAnalysis (pp-only): kn=${kn}, e_com=${e_com}"

root -b -l << EOF
	.L PythiaAnalysisClasses.h
	PythiaTruthAnalysis pw($kn, (bool)0, $e_com, (bool)0, (bool)1);
	pw.Run();
	.q;
EOF
