#!/bin/bash

# The argument gets passed as the argument $(Process) in batch-submission script
file_batch=$(( $1 + 1 ))

cd $PWD

# Setup ATLAS environment
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

# Setup LCG
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

# Run the analysis
root -b -l << EOF
	.L PPDataNTupleFirstPass.c

	PPDataNTupleFirstPass pp_17_MB (17, $file_batch);
	pp_17_MB.isMinBias = true;
	// pp_17_MB.resonance_cut_mode = 0 // uncomment for no-resn-cut
	pp_17_MB.Run();

	.q;

EOF

