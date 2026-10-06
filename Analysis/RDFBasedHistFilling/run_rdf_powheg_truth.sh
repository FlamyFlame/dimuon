#!/usr/bin/env bash
set -eo pipefail

set +e
set +u
source ~/setup.sh
setup_status=$?
set -e
set -u

if [[ ${setup_status} -ne 0 ]]; then
  echo "Environment setup failed (source ~/setup.sh)."
  exit ${setup_status}
fi

# Optional arg: "nores" -> bb no-resonance-cut NTuple only (signal acceptance), else nominal.
WITH_RES_CUTS=true
if [[ "${1:-}" == "nores" ]]; then WITH_RES_CUTS=false; fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

root -l -b <<EOF
.L RDFBasedHistFillingPowhegTruth.cxx+
RunRDFBasedHistFillingPowhegTruth(${WITH_RES_CUTS});
.q
EOF
