#!/bin/bash
file_batch=$(( $1 + 1 ))
cd $PWD
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"
root -b -l << ROOTEOF
	.L DataAnalysisClasses.h
	PPAnalysis pp_24 (24, $file_batch);
	pp_24.trigger_mode = 3;
	pp_24.resonance_cut_mode = 0;
	pp_24.Run();
	.q;
ROOTEOF
