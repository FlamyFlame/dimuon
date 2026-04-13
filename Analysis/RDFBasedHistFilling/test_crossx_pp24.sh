#!/usr/bin/env bash
# Test FillHistogramsCrossx for pp24 with trigger_mode=3 (2mu4)

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

echo "[TEST] Running FillHistogramsCrossx for pp24 with trigger_mode=3 (2mu4)..."

root -l -b <<EOF
// Load and compile PP (+ forces recompilation)
.L RDFBasedHistFillingPP.cxx+

// Create PP instance for run_year=24
RDFBasedHistFillingPP pp(24);

// Configure for trigger_mode=3 (2mu4, crossx measurements for pp)
pp.trigger_mode = 3;

// Set mindR_trig to -1 to disable mindR filtering
pp.mindR_trig   = -1;

std::cout << "\\n[TEST] Starting pp24 FillHistogramsCrossx run..." << std::endl;
pp.Run();
std::cout << "\\n[TEST] pp24 FillHistogramsCrossx completed successfully!" << std::endl;

.q
EOF

echo "[TEST] pp24 test completed. Output file:"
ls -lh /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024*.root 2>/dev/null | tail -1
