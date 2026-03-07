#!/usr/bin/env bash
# Run RDF-based histogram filling for PbPb 2023 legacy data.
#
# "Legacy" means the skimmed NTuple was produced WITHOUT the new dimuon trigger
# matching features added in the 2025 skimming update:
#   - No per-muon leg information for the asymmetric HLT_mu4_mu4noL1_L1MU3V trigger
#   - No minimum-dR (mindR) requirement applied during skimming
#
# Consequently:
#   useMu4NoL1Leg = false  -- do not require the 2nd muon to pass the mu4noL1 leg
#   mindR_trig    = -1     -- search for input files WITHOUT a _mindR_X_XX suffix
#
# Load order: RDFBasedHistFillingPowheg.cxx is loaded first (without recompiling;
# ROOT picks up the existing RDFBasedHistFillingPowheg_cxx.so).  This pre-populates
# the base-class template instantiations that cling needs to JIT-link PbPb.
# RDFBasedHistFillingBaseClass.cxx has #pragma once, so it is not redefined.
#
# Usage:
#   ./run_pbpb23_hist_filling_legacy.sh [include_upc|default]
#
# Defaults:
#   ctr_binning = include_upc

set -eo pipefail

ctr_binning="${1:-include_upc}"

case "${ctr_binning}" in
  include_upc|default) ;;
  *)
    echo "Invalid ctr_binning: ${ctr_binning}. Use include_upc|default"
    exit 1
    ;;
esac

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
// Pre-load Powheg .so so that RDFBasedHistFillingBaseClass symbols are already
// available; PbPb.cxx then skips the base class via #pragma once.
.L RDFBasedHistFillingPowheg.cxx

// Load PbPb without ACLiC compilation (no +).
.L RDFBasedHistFillingPbPb.cxx

RDFBasedHistFillingPbPb pbpb(23, "${ctr_binning}");

// Legacy data: no mu4noL1-leg info, no mindR-filtered skimming -> use old file names
pbpb.useMu4NoL1Leg = false;  // do not filter on 2nd muon passing the mu4noL1 leg
pbpb.mindR_trig    = -1;     // negative -> search for input files without _mindR_X_XX suffix

pbpb.trigger_mode  = 1;      // single-mu4 triggered data (dimuon trigger efficiency mode)

pbpb.Run();
.q
EOF
