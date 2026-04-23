#!/bin/bash

file_batch=$(( $1 + 1 ))

cd $PWD

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

root -b -l << EOF
	.L DataAnalysisClasses.h

	PbPbAnalysis pbpb_24 (24, $file_batch);
	pbpb_24.pbpb24_mu4_NO_trig_calc = true;
	pbpb_24.nevents_max = 1000;
	pbpb_24.is_test_run = true;
	pbpb_24.Run();

	.q;

EOF
