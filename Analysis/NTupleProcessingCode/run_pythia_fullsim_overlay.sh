#!/bin/bash

# Run Pythia fullsim HIJING overlay dimuon analysis.
# Processes all available kn NTUP files in pythia_fullsim_hijing_overlay_test_sample.
# Usage: ./run_pythia_fullsim_overlay.sh

cd $PWD

# Setup ATLAS environment
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

# Setup LCG
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

# Run the analysis
root -b -l << EOF
	.L PythiaAnalysisClasses.h

	PythiaFullSimOverlayAnalysis py(FullSimSampleType::hijing);
	// TEST sample: reads the overlay test sample AND uses the pp beam only (only beam produced).
	// ONE switch -- it also selects the input dir. See FullSimSampleType.h.
	py.isTestSample = true;
	// AMI PROVENANCE: the overlay is built on the pp-beam evgen DSIDs. AMI files are keyed by
	// beam+slice only, so declare them; InitInputFullsim throws if another production is read.
	py.expected_ami_dsids = {802776, 802777, 802778, 802779, 802780, 802781};
	// py.nevents_max = 1000;  // uncomment to test with limited events
	// py.debug_mode = true;   // uncomment for verbose debug output
	py.fill_kn_trees_fullsim = true;
	py.Run();

	.q;

EOF
