#!/bin/bash

# Setup ATLAS environment
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

# Setup LCG
lsetup "views LCG_105 x86_64-centos7-gcc11-opt"

# Run the analysis

root -b -l << EOF
	.L PbPbDataNTupleFirstPass.c

	PbPbDataNTupleFirstPass pbpb_run2;
	pbpb_run2.isRun3 = false;
	pbpb_run2.Run();

	.q;
EOF

