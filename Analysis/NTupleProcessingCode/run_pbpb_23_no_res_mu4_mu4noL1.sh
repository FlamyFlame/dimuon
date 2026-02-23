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
	.L DataAnalysisClasses.h

	PbPbAnalysis pbpb_23 (23, $file_batch);
 	pbpb_23.resonance_cut_mode = 0;
 	pbpb_23.trigger_mode = 2;
	pbpb_23.Run();

	.q;

EOF

