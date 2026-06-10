#!/bin/bash

# Run Pythia pp24 fullsim: single-muon tree output for reco efficiency sanity check.
# Usage: ./run_pythia_fullsim_single_muon.sh

cd "$(dirname "$0")"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

root -b -l << EOF
	.L PythiaAnalysisClasses.h

	PythiaFullSimAnalysis py(0);
	py.fill_kn_trees_fullsim = true;
	py.output_single_muon_tree = true;
	py.extra_output_suffix = "_single_muon";
	py.Run();

	.q;

EOF
