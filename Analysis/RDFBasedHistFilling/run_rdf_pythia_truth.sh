#!/usr/bin/env bash
set -eo pipefail

# Usage:
#   ./run_rdf_pythia_truth.sh [private|nonprivate] [5.36|5.02] [withcuts|nocuts]
#
# Defaults:
#   mode    = nonprivate
#   ecom    = 5.36
#   cuts    = nocuts

mode="${1:-nonprivate}"
ecom="${2:-5.36}"
cuts="${3:-nocuts}"

is_private=false
case "${mode}" in
  private)    is_private=true ;;
  nonprivate) is_private=false ;;
  *)
    echo "Invalid mode: ${mode}. Use private|nonprivate"
    exit 1
    ;;
esac

with_cuts=false
case "${cuts}" in
  withcuts) with_cuts=true ;;
  nocuts)   with_cuts=false ;;
  *)
    echo "Invalid cuts option: ${cuts}. Use withcuts|nocuts"
    exit 1
    ;;
esac

if [[ "${ecom}" != "5.36" && "${ecom}" != "5.02" ]]; then
  echo "Invalid E_COM: ${ecom}. Use 5.36 or 5.02"
  exit 1
fi

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
.L RDFBasedHistFillingPythiaTruth.cxx+
RunRDFBasedHistFillingPythiaTruth(${is_private}, ${ecom}, ${with_cuts});
.q
EOF
