#!/usr/bin/env bash
set -euo pipefail

set +u
source ~/setup.sh
set -u

ANALYSIS_DIR="/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
PLOT_DIR="${ANALYSIS_DIR}/plotting_codes/single_b_analysis"
OUT_DIR="/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/single_b_analysis"

mkdir -p "${OUT_DIR}"

run_and_log() {
  local label="$1"
  local cmd="$2"
  echo "[RUN] ${label}"
  bash -lc "${cmd}"
  echo "[OK ] ${label}"
}

# 1) Fill crossx histograms
# PbPb runs first: test_crossx_pp24.sh loads the PbPb .so in a separate session to resolve
# the dynamic_cast<RDFBasedHistFillingPbPb*> typeinfo dependency before linking PP.
run_and_log "RDF crossx pbpb23" "cd '${RDF_DIR}' && bash test_crossx_pbpb23.sh"
run_and_log "RDF crossx pbpb24" "cd '${RDF_DIR}' && bash test_crossx_pbpb24.sh"
run_and_log "RDF crossx pp24" "cd '${RDF_DIR}' && bash test_crossx_pp24.sh"

# 2) Plot
run_and_log "Plot crossx pp24" "cd '${ANALYSIS_DIR}' && root -l -b -q 'plotting_codes/single_b_analysis/plot_single_b_crossx_pp.cxx'"
run_and_log "Plot crossx pbpb combined" "cd '${ANALYSIS_DIR}' && root -l -b -q 'plotting_codes/single_b_analysis/plot_single_b_crossx_pbpb.cxx'"

# 3) Validate png outputs are non-empty
validate_png() {
  local png="$1"
  if [[ ! -f "${png}" ]]; then
    echo "[ERR] Missing PNG: ${png}" >&2
    exit 1
  fi
  if [[ ! -s "${png}" ]]; then
    echo "[ERR] Empty PNG: ${png}" >&2
    exit 1
  fi
  echo "[OK ] ${png} ($(stat -c%s "${png}") bytes)"
}

validate_png "${OUT_DIR}/pp24/pp24_crossx_pair_pt_pair_eta.png"
validate_png "${OUT_DIR}/pp24/pp24_crossx_pair_pt_minv.png"
validate_png "${OUT_DIR}/pp24/pp24_crossx_pair_pt_dr.png"
validate_png "${OUT_DIR}/pp24/pp24_crossx_pair_pt_in_eta_subplots_dr_lines.png"

# Combined PbPb (discovers which years are available; check a representative subset)
COMB_TAA="${OUT_DIR}/pbpb_23_24_combined/TAA_weighted"
for ctr in ctr0_5 ctr5_10 ctr10_20 ctr20_30 ctr30_50 ctr50_80; do
  validate_png "${COMB_TAA}/pbpb_combined_${ctr}_pair_pt_pair_eta.png"
  validate_png "${COMB_TAA}/pbpb_combined_${ctr}_pair_pt_in_eta_subplots_dr_lines.png"
done

echo "[DONE] All crossx filling and plotting finished successfully."
