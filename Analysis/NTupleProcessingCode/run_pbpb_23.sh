#!/bin/bash

# The argument gets passed as the argument $(Process) in batch-submission script
file_batch=$(( $1 + 1 ))

cd $PWD

# Setup ATLAS environment
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

# Setup LCG
lsetup "views LCG_105 x86_64-centos7-gcc11-opt"

# Run the analysis
root -b -l << EOF
	.L PbPbDataNTupleFirstPass.c

	PbPbDataNTupleFirstPass pbpb_23;
	pbpb_23.run_year = 23;
	pbpb_23.file_batch = $file_batch;
	pbpb_23.Run();

	.q;

EOF

