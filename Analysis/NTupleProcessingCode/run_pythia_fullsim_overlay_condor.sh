#!/bin/bash
# Condor worker: Pythia fullsim overlay NTuple processing.
# Usage: ./run_pythia_fullsim_overlay_condor.sh <sample_type> [overlay_year]
#   sample_type : hijing | zmumu | data
#   overlay_year: 24 (default, pbpb24 test sample) | 23 (_pbpb23 test sample); hijing only
#
# TEMPORARY: Currently processes all 6 kn ranges in a single job
# (fullsim test sample is small enough). When full samples become
# available, this script and the .sub file will need to be updated
# to support batched kn processing and multiple Condor jobs.

sample_type="${1:-hijing}"
# HIJING-overlay conditions year (2nd argument, default 24; 23 = the _pbpb23 sample) -- see
# FullSimSampleType.h. Drives the input directory AND the output label together.
overlay_year="${2:-24}"

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
    .L PythiaAnalysisClasses.h+

    PythiaFullSimOverlayAnalysis py(FullSimSampleType::${sample_type});
    // TEST sample: reads the overlay test sample AND uses the pp beam only (only beam produced).
    // ONE switch -- it also selects the input dir. See FullSimSampleType.h.
    py.isTestSample = true;
    py.overlay_pbpb_year = ${overlay_year};
    // AMI PROVENANCE: the overlay is built on the pp-beam evgen DSIDs. AMI files are keyed by
    // beam+slice only, so declare them; InitInputFullsim throws if another production is read.
    py.expected_ami_dsids = {802776, 802777, 802778, 802779, 802780, 802781};
    py.fill_kn_trees_fullsim = true;
    py.Run();

    .q;

EOF
