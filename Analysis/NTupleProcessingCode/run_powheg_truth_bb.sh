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
	.L PowhegAnalysisClasses.h 

	PowhegTruthAnalysis pw($file_batch, "bb");
	pw.Run();

	.q;

EOF

