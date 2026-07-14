#!/bin/bash

# Run Pythia fullsim HIJING overlay with trigger propagation (store_mc_trigger):
# pair trees with m1/m2.passmu4 + pass2mu4 from the trigger-enabled _July2026 skims.
# Output suffix: _mc_trig (nominal outputs untouched).
# Usage: ./run_pythia_fullsim_overlay_mc_trig.sh

cd "$(dirname "$0")"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

root -b -l << EOF
	.L PythiaAnalysisClasses.h+

	PythiaFullSimOverlayAnalysis py(FullSimSampleType::hijing);
	// TEST sample: reads the overlay test sample AND uses the pp beam only (only beam produced).
	// ONE switch -- it also selects the input dir. See FullSimSampleType.h.
	py.isTestSample = true;
	// AMI PROVENANCE: the overlay is built on the pp-beam evgen DSIDs. AMI files are keyed by
	// beam+slice only, so declare them; InitInputFullsim throws if another production is read.
	py.expected_ami_dsids = {802776, 802777, 802778, 802779, 802780, 802781};
	py.fill_kn_trees_fullsim = true;
	py.store_mc_trigger = true;
	py.Run();

	.q;

EOF
