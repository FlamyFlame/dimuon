#!/usr/bin/env bash
# Low-mass dimuon TEMPLATE-FIT data pass for pp24 (trigger_mode=3, 2mu4).
# Reads the _no_res_cut ntuples (resonances PRESENT) and writes the 0-4 GeV minv
# template spectra D_OS/D_SS (1D + 2D vs pair pT/eta) to a DISTINCT output
# (histograms_real_pairs_pp_2024_2mu4_nominal_template_fit.root). Part of the crossx
# pipeline: feeds the low-mass template fit -> R_AA (low_mass_dimuon_template_fit.md).

set -Eeuo pipefail
set +e; set +u
source ~/setup.sh
setup_status=$?
set -e; set -u
if [[ ${setup_status} -ne 0 ]]; then echo "Environment setup failed (source ~/setup.sh)."; exit ${setup_status}; fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

echo "[RUN] pp24 low-mass template-fit pass (trigger_mode=3, _no_res_cut)..."

root -l -b <<EOF
.L RDFBasedHistFillingPP.cxx+
RDFBasedHistFillingPP pp(24);
pp.trigger_mode = 3;              // 2mu4 (same as nominal crossx)
pp.low_mass_template_calc = true; // -> _no_res_cut input + _template_fit output; minv-only fill
pp.output_generic_hists = false;  // light pass: weights computed inline in the template block
pp.output_gapcut_hists  = false;
std::cout << "\\n[RUN] Starting pp24 template-fit run..." << std::endl;
pp.Run();
std::cout << "\\n[RUN] pp24 template-fit completed successfully!" << std::endl;
gSystem->Exit(0);
EOF

echo "[RUN] pp24 template-fit completed. Output file:"
ls -lh /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024*template_fit*.root 2>/dev/null | tail -1
