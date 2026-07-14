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
	// Reads the TEST sample; pp_only=true forces the pp beam alone (weight 1) as a
	// cross-check against that sample's nominal 4-beam isospin average.
	py.isTestSample = true;
	// AMI PROVENANCE: declare the TEST sample's pp-beam DSIDs (see PythiaAlgCoreT.h).
	py.expected_ami_dsids = {802776, 802777, 802778, 802779, 802780, 802781};
	py.fill_kn_trees_fullsim = true;
	py.Run();
	.q;
EOF
