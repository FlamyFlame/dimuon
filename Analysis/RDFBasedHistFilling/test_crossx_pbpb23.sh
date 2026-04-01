#!/usr/bin/env bash
# Test FillHistogramsCrossx for pbpb23 with trigger_mode=2

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

echo "[TEST] Running FillHistogramsCrossx for pbpb23 with trigger_mode=2..."

root -l -b <<EOF
// Pre-load Powheg .so for base class symbols
.L RDFBasedHistFillingPowheg.cxx

// Load and compile PbPb (+ forces recompilation)
.L RDFBasedHistFillingPbPb.cxx+

// Create PbPb instance for run_year=23
RDFBasedHistFillingPbPb pbpb(23);

// Configure for trigger_mode=2 (crossx measurements)
pbpb.trigger_mode = 2;

// Set mindR_trig to -1 to disable mindR filtering
pbpb.mindR_trig   = -1;

std::cout << "\\n[TEST] Starting pbpb23 FillHistogramsCrossx run..." << std::endl;
pbpb.Run();
std::cout << "\\n[TEST] pbpb23 FillHistogramsCrossx completed successfully!" << std::endl;

.q
EOF

echo "[TEST] pbpb23 test completed. Output file:"
ls -lh /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/histograms_real_pairs_pbpb_2023*.root 2>/dev/null | tail -1
