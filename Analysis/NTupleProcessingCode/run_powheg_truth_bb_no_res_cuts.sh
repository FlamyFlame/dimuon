#!/bin/bash

# The argument gets passed as the argument $(Process) in batch-submission script
file_batch=$(( $1 + 1 ))
use_local=${2:-0}

cd $PWD

# Setup ATLAS environment
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

# Setup LCG
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

# Run the analysis -- NO-RESONANCE-CUT mode (truth signal acceptance needs ALL masses;
# docs/tracking/signal_acceptance_mass_window.md R2). Otherwise identical to run_powheg_truth_bb.sh.
root -b -l << EOF
	.L PowhegAnalysisClasses.h 

	PowhegTruthAnalysis pw($file_batch, "bb", $use_local);
	pw.setTurnDataResonanceCutsOn(false);   // keep resonance-tagged pairs -> *_no_data_resonance_cuts_part<N>.root
	pw.Run();

	.q;

EOF

