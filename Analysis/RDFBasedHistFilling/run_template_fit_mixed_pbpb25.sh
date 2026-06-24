#!/usr/bin/env bash
# Low-mass dimuon TEMPLATE-FIT data pass for pbpb25 (trigger_mode=1, single_mu4, nominal mode).
# Reads _no_res_cut (resonances PRESENT) -> per-centrality 0-4 GeV minv template spectra
# D_OS/D_SS (1D + 2D vs pair pT/eta) to DISTINCT output
# (histograms_real_pairs_pbpb_2025_single_mu4_no_trg_plots_nominal_template_fit.root).
set -Eeuo pipefail
set +e; set +u
source ~/setup.sh
setup_status=$?
set -e; set -u
if [[ ${setup_status} -ne 0 ]]; then echo "Environment setup failed."; exit ${setup_status}; fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
echo "[RUN] pbpb25 mixed-event template (T_mix) pass (trigger_mode=1, _scrambled, mixed-event)..."
root -l -b <<ROOTEOF
.L RDFBasedHistFillingPbPb.cxx+
RDFBasedHistFillingPbPb pbpb(25);
pbpb.trigger_mode = 1;                       // single_mu4 (same as nominal crossx)
pbpb.mu4_nominal_pbpb_NO_trig_calc = true;   // nominal muon selection (no trig-eff)
pbpb.mindR_trig = 0.02;
pbpb.low_mass_template_calc = true;
  pbpb.mixed_event_template = true;          // -> _scrambled input + _template_fit_mixed_event output
pbpb.output_generic_hists = false;           // light pass: weights inline in the template block
std::cout << "\\n[RUN] Starting pbpb25 template-fit run..." << std::endl;
pbpb.Run();
std::cout << "\\n[RUN] pbpb25 mixed-event template completed successfully!" << std::endl;
.q
ROOTEOF
echo "[RUN] pbpb25 mixed-event template completed. Output file:"
ls -lh /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2025/histograms_real_pairs_pbpb_2025*template_fit*.root 2>/dev/null | tail -1
