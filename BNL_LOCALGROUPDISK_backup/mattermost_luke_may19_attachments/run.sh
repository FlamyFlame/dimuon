#!/bin/bash
set -euo pipefail

WILD="$1"
ARG1="$1"
ARG2="$2"
ARG3="$3"
# `Process` (the job number) is propagated to `JOBNO``
    # (You might want `JOBNO` to propagate as an argument to your macro, or not)
JOBNO=$4

# Make the output file directory if it doesn't exist
mkdir -p "$PWD/output_files"

# For convenience define OUTDIR as the output file directory
OUTDIR="output_files"

# Export the user proxy
export X509_USER_PROXY="$PWD/x509up_u102136"

# Setup root on lxplus (change the setup script as necessary)
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_108_ATLAS_2/x86_64-el9-gcc14-opt/setup.sh
set -u

# Run your macro
root -b -l <<EOF
.L myMacro.C
myMacro("${WILD}", "${ARG1}", "${ARG2}", "${ARG3}")
.q
EOF