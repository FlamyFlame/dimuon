#!/bin/bash

# Run Pythia pp24 fullsim dimuon analysis, pp isospin only.
# Output: muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_pp_only.root
# Usage: ./run_pythia_fullsim_pp_only.sh

cd "$(dirname "$0")"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

root -b -l << EOF
	.L PythiaAnalysisClasses.h
	PythiaFullSimAnalysis py(0, (bool)0, (bool)1);
	py.fill_kn_trees_fullsim = true;
	py.Run();
	.q;
EOF
