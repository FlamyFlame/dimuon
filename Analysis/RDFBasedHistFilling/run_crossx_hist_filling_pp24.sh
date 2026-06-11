#!/usr/bin/env bash
# Fill crossx histograms for pp24 with trigger_mode=3 (2mu4)

set -Eeuo pipefail

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

echo "[RUN] Running FillHistogramsCrossx for pp24 with trigger_mode=3 (2mu4)..."

root -l -b <<EOF
// Load and compile PP
.L RDFBasedHistFillingPP.cxx+

// Create PP instance for run_year=24
RDFBasedHistFillingPP pp(24);

// Configure for trigger_mode=3 (2mu4, crossx measurements for PP24)
pp.trigger_mode = 3;
pp.output_generic_hists = true;
pp.output_gapcut_hists = true;

std::cout << "\\n[RUN] Starting pp24 FillHistogramsCrossx run..." << std::endl;
pp.Run();
std::cout << "\\n[RUN] pp24 FillHistogramsCrossx completed successfully!" << std::endl;

gSystem->Exit(0);
EOF

echo "[RUN] pp24 hist filling completed. Output file:"
ls -lh /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024*.root 2>/dev/null | tail -1
