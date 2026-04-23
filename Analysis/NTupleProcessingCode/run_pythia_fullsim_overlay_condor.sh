#!/bin/bash
# Condor worker: Pythia fullsim overlay NTuple processing.
# Usage: ./run_pythia_fullsim_overlay_condor.sh <sample_type>
#   sample_type: hijing | zmumu | data
#
# TEMPORARY: Currently processes all 6 kn ranges in a single job
# (fullsim test sample is small enough). When full samples become
# available, this script and the .sub file will need to be updated
# to support batched kn processing and multiple Condor jobs.

sample_type="${1:-hijing}"

case "$sample_type" in
    hijing|zmumu|data) ;;
    *) echo "Unknown sample_type: $sample_type (use hijing|zmumu|data)"; exit 1 ;;
esac

cd "$(dirname "$0")"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

echo "Starting PythiaFullSimOverlayAnalysis: sample_type=${sample_type}"

root -b -l << EOF
    .L PythiaAnalysisClasses.h

    PythiaFullSimOverlayAnalysis py(FullSimSampleType::${sample_type});
    py.fill_kn_trees_fullsim = true;
    py.Run();

    .q;

EOF
