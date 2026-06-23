#!/bin/bash
file_batch=$(( $1 + 1 ))
cd $PWD
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"
root -b -l << ROOTEOF
	.L DataAnalysisClasses.h
	PbPbAnalysis pbpb_24 (24, $file_batch);
	pbpb_24.pbpb_run3_mu4_force_nominal = true;
	pbpb_24.resonance_cut_mode = 0;
	pbpb_24.Run();
	.q;
ROOTEOF
