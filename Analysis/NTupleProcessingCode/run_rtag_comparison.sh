#!/bin/bash
# Run Pythia fullsim overlay NTuple processing for a single r-tag variant.
# Usage: ./run_rtag_comparison.sh <suffix> <input_dir> [pythia_only_bc]
#   suffix:    _rtag_orig | _rtag_signalOnlyTruth | _rtag_no_overlay
#   input_dir: path to directory containing the symlinked NTUP
#   pythia_only_bc: if "1", restrict barcode cache to Pythia-only truth

SUFFIX="$1"
INPUT_DIR="$2"
PYTHIA_ONLY_BC="${3:-0}"

if [[ -z "$SUFFIX" || -z "$INPUT_DIR" ]]; then
    echo "Usage: $0 <suffix> <input_dir> [pythia_only_bc]"
    exit 1
fi

cd "$PWD"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

root -b -l << EOF
    .L PythiaAnalysisClasses.h

    PythiaFullSimOverlayAnalysis py(FullSimSampleType::hijing);
    py.fullsim_input_dir_override = "${INPUT_DIR}";
    py.extra_output_suffix = "${SUFFIX}";
    py.fill_kn_trees_fullsim = true;
    py.pythia_only_barcode_cache = ${PYTHIA_ONLY_BC};
    py.Run();

    .q;

EOF
