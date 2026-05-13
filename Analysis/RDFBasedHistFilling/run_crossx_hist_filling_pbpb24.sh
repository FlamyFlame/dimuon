#!/usr/bin/env bash
# Fill crossx histograms for pbpb24 with trigger_mode=1 (mu4)

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

echo "[RUN] Running FillHistogramsCrossx for pbpb24 with trigger_mode=1 (mu4)..."

root -l -b <<EOF
// Load and compile PbPb (+ forces recompilation)
.L RDFBasedHistFillingPbPb.cxx+

// Create PbPb instance for run_year=24
RDFBasedHistFillingPbPb pbpb(24);

// Configure for trigger_mode=1 (mu4 single muon crossx)
// doTrigEffcy left at default (true); run_crossx=true prevents trig effcy filling
// so the output filename stays _single_mu4 (no _no_trg_plots suffix)
pbpb.trigger_mode = 1;

// Use mindR_trig=0.02 to read latest ntuples with event selection
pbpb.mindR_trig   = 0.02;

std::cout << "\\n[RUN] Starting pbpb24 FillHistogramsCrossx run..." << std::endl;
pbpb.Run();
std::cout << "\\n[RUN] pbpb24 FillHistogramsCrossx completed successfully!" << std::endl;

.q
EOF

echo "[RUN] pbpb24 hist filling completed. Output file:"
ls -lh /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2024/histograms_real_pairs_pbpb_2024*.root 2>/dev/null | tail -1
