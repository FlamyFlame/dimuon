#!/bin/bash

# Run Pythia pp24 fullsim with trigger propagation (store_mc_trigger):
# pair trees with m1/m2.passmu4 + pass2mu4 from the trigger-enabled _July2026 skims.
# Output suffix: _mc_trig (nominal outputs untouched). Files without trigger
# branches (old skims) are skipped with a loud message.
# Usage: ./run_pythia_fullsim_mc_trig.sh

cd "$(dirname "$0")"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

root -b -l << EOF
	.L PythiaAnalysisClasses.h+

	PythiaFullSimAnalysis py(0);
	py.fill_kn_trees_fullsim = true;
	py.store_mc_trigger = true;
	py.Run();

	.q;

EOF
