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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

root -l -b <<EOF
.L RDFBasedHistFillingPowhegTruth.cxx+
RunRDFBasedHistFillingPowhegTruth();
.q
EOF
