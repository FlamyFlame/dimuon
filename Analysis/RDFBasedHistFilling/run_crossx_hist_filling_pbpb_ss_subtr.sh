#!/usr/bin/env bash
# OS-SS subtraction crossx hist filling for ONE Pb+Pb year (mu4, trigger_mode=1), usage:
#   ./run_crossx_hist_filling_pbpb_ss_subtr.sh <yy>      # yy in 23 24 25 26
# Writes histograms_real_pairs_pbpb_20<yy>_<trig>_no_trg_plots_nominal_ss_subtr_crossx.root (OS families under
# their nominal names + SS twins). The NOMINAL crossx file is never touched. Env: BUILD_DIR = ACLiC build dir
# (default: a private one, so a concurrent session's .so is not rebuilt underneath it).
set -eo pipefail
yr="${1:?usage: $0 <yy>}"
set +e; set +u; source ~/setup.sh; setup_status=$?; set -e; set -u
[[ ${setup_status} -eq 0 ]] || { echo "Environment setup failed (source ~/setup.sh)."; exit ${setup_status}; }
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-/tmp/ss_subtr_build_${yr}}"; mkdir -p "${BUILD_DIR}"
stamp="$(mktemp)"
echo "[RUN] FillHistogramsCrossx OS-SS subtraction pass, pbpb${yr} (mu4), build dir ${BUILD_DIR}"
root -l -b <<ROOTEOF
gSystem->SetBuildDir("${BUILD_DIR}", true);
.L RDFBasedHistFillingPbPb.cxx+
RDFBasedHistFillingPbPb pbpb(${yr});
pbpb.trigger_mode = 1;
pbpb.mu4_nominal_pbpb_NO_trig_calc = true;
pbpb.mindR_trig   = 0.02;
pbpb.ss_subtr_crossx = true;
pbpb.Run();
std::cout << "\n[RUN] pbpb${yr} OS-SS subtraction fill completed" << std::endl;
.q
ROOTEOF
out=$(ls /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_20${yr}/histograms_real_pairs_pbpb_20${yr}*_ss_subtr_crossx.root | head -1)
# ROOT swallows exceptions thrown inside the RDF loop and exits 0: require a FRESH file with FILLED OS and SS twins.
[[ "$out" -nt "$stamp" ]] || { echo "FAIL: $out is older than this run (the fill threw)"; exit 1; }
root -l -b >/dev/null 2>&1 <<CHK || { echo "FAIL: $out lacks a filled OS/SS pair"; exit 1; }
TFile* f = TFile::Open("${out}", "READ");
if (!f || f->IsZombie()) gSystem->Exit(2);
for (const char* n : {"h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_ctr0_5", "h2d_ss_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_ctr0_5"}) {
  TH1* h = (TH1*)f->Get(n);
  if (!h || h->GetEntries() <= 0) gSystem->Exit(3);
}
gSystem->Exit(0);
CHK
rm -f "$stamp"
ls -lh "$out"
