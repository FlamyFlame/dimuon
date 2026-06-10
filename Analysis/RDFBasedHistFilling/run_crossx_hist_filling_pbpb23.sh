#!/usr/bin/env bash
# Fill crossx histograms for pbpb23 with trigger_mode=1 (mu4)

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

echo "[RUN] Running FillHistogramsCrossx for pbpb23 with trigger_mode=1 (mu4)..."

root -l -b <<EOF
// Load and compile PbPb (+ forces recompilation)
.L RDFBasedHistFillingPbPb.cxx+

// Create PbPb instance for run_year=23
RDFBasedHistFillingPbPb pbpb(23);

// Configure for trigger_mode=1 (single_mu4, crossx/nominal measurements for PbPb23)
pbpb.trigger_mode = 1;
pbpb.mu4_nominal_pbpb_NO_trig_calc = true;
pbpb.useCoarseQEtaBin = true;

// Use mindR_trig=0.02 to read latest ntuples with event selection
pbpb.mindR_trig   = 0.02;

std::cout << "\\n[RUN] Starting pbpb23 FillHistogramsCrossx run..." << std::endl;
pbpb.Run();
std::cout << "\\n[RUN] pbpb23 FillHistogramsCrossx completed successfully!" << std::endl;

.q
EOF

echo "[RUN] pbpb23 hist filling completed. Output file:"
ls -lh /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/histograms_real_pairs_pbpb_2023*.root 2>/dev/null | tail -1
