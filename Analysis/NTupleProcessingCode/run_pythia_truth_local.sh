#!/bin/bash
# Run Pythia truth analysis in local-file mode for the 5.02 TeV 40-70 GeV bin.
# Usage: ./run_pythia_truth_local.sh [kn_index]
#   kn_index: must be 3 (40-70 GeV), default 3

kn=${1:-3}
is_private=0
e_com=5.02
use_local=1

if [[ "$kn" -ne 3 ]]; then
    echo "ERROR: local mode only supports kn=3 (40-70 GeV)."
    echo "Use pnfs mode for other kn bins (run_pythia_truth_kn.sh <kn> 0 5.02 0)."
    exit 1
fi

batch_num=$kn

cd "$(dirname "$0")"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

echo "Starting local PythiaTruthAnalysis: kn=${kn}, batch_num=${batch_num}, e_com=${e_com}, use_local=${use_local}"

root -b -l << EOF
    .L PythiaAnalysisClasses.h
    PythiaTruthAnalysis pw($batch_num, (bool)$is_private, $e_com, (bool)$use_local);
    pw.nevents_max = 1000; // Optional quick test mode
    pw.Run();
    .q;
EOF
