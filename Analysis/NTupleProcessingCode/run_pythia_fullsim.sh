#!/bin/bash

# Run Pythia pp24 fullsim dimuon analysis.
# Processes all available kn×beam NTUP files in pythia_fullsim_test_sample.
# Usage: ./run_pythia_fullsim.sh

cd $PWD

# Setup ATLAS environment
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

# Setup LCG
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

# Run the analysis
root -b -l << EOF
	.L PythiaAnalysisClasses.h

	PythiaFullSimAnalysis py(0);
	// TEST sample: reads pythia_fullsim_test_sample/ AND uses the 4-isospin-beam average
	// (that sample was produced with 4 beams by mistake). ONE switch -- see FullSimSampleType.h.
	py.isTestSample = true;
	// AMI PROVENANCE: AMI files are keyed by beam+slice only, so nothing otherwise stops
	// this run from being weighted with ANOTHER production's cross-sections. Declare the
	// TEST sample's DSIDs; InitInputFullsim throws on a mismatch.
	py.expected_ami_dsids = {802758, 802759, 802760, 802761, 802762, 802763,
                                802764, 802765, 802766, 802767, 802768, 802769,
                                802770, 802771, 802772, 802773, 802774, 802775,
                                802776, 802777, 802778, 802779, 802780, 802781};
	// py.nevents_max = 1000;  // uncomment to test with limited events
	// py.debug_mode = true;   // uncomment for verbose debug output
	py.fill_kn_trees_fullsim = true;
	py.Run();

	.q;

EOF
