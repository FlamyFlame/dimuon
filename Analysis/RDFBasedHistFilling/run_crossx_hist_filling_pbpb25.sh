#!/usr/bin/env bash
# Fill crossx histograms for pbpb25 with trigger_mode=1 (mu4)

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

echo "[RUN] Running FillHistogramsCrossx for pbpb25 with trigger_mode=1 (mu4)..."

root -l -b <<EOF
// Load and compile PbPb (+ forces recompilation)
.L RDFBasedHistFillingPbPb.cxx+

// Create PbPb instance for run_year=25
RDFBasedHistFillingPbPb pbpb(25);

// Configure for trigger_mode=1 (single_mu4, crossx measurements for PbPb25)
pbpb.trigger_mode = 1;

// Use mindR_trig=0.02 to read latest ntuples with event selection
pbpb.mindR_trig   = 0.02;

std::cout << "\\n[RUN] Starting pbpb25 FillHistogramsCrossx run..." << std::endl;
pbpb.Run();
std::cout << "\\n[RUN] pbpb25 FillHistogramsCrossx completed successfully!" << std::endl;

.q
EOF

echo "[RUN] pbpb25 hist filling completed. Output file:"
ls -lh /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2025/histograms_real_pairs_pbpb_2025*.root 2>/dev/null | tail -1
