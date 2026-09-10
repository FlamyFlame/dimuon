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

# ─── PREREQUISITE ──────────────────────────────────────────────────────────────────────────────
# This driver fills and plots the cross-sections ONLY. It does NOT run the trigger-efficiency
# stage, and the Pb+Pb crossx fill CONSUMES the single-muon turn-on TF1 fits
# (<year>/trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit.root). If those are from a
# previous round, the fill either throws inside the RDF event loop -- where ROOT swallows the
# exception and still exits 0 -- or silently weights every pair with the wrong efficiency.
# Run pipelines/run_pbpb_all.sh (trig-eff -> medium -> crossx, serialized) for a full refresh;
# use this driver only when the fits on disk are known to be current.
# ───────────────────────────────────────────────────────────────────────────────────────────────

# Everything this run validates must be NEWER than this marker. Existence + non-zero size passes
# on every leftover from the previous round, and a plotting macro that throws still exits 0, so
# without a freshness test the script can print "All crossx filling and plotting finished
# successfully" having produced nothing at all.
RUN_MARKER="$(mktemp)"
trap 'rm -f "${RUN_MARKER}"' EXIT

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
# Every Pb+Pb data-taking year, filtered on input presence: a year whose muon_pairs
# ntuple is not on disk yet (its skim or NTuple processing still running) is announced and
# skipped, rather than failing the run for the years that ARE ready.  PBPB_YEARS is also
# what names the combined output directory below, so the two cannot drift apart.
DATA_BASE="/usatlas/u/yuhanguo/usatlasdata/dimuon_data"
PBPB_YEARS=()
for yr in 23 24 25 26; do
  if compgen -G "${DATA_BASE}/pbpb_20${yr}/muon_pairs_pbpb_20${yr}_single_mu4"*.root > /dev/null; then
    PBPB_YEARS+=("${yr}")
  else
    echo "[SKIP] Pb+Pb 20${yr}: no muon_pairs ntuple on disk -- crossx filling skipped." >&2
  fi
done
if [[ ${#PBPB_YEARS[@]} -eq 0 ]]; then
  echo "[FATAL] no Pb+Pb year has a muon_pairs ntuple on disk." >&2
  exit 1
fi
echo "[INFO] Pb+Pb crossx years: ${PBPB_YEARS[*]}"
for yr in "${PBPB_YEARS[@]}"; do
  run_and_log "RDF crossx pbpb${yr}" "cd '${RDF_DIR}' && bash run_crossx_hist_filling_pbpb${yr}.sh"
done
# NOT test_crossx_pp24.sh. Both write the SAME nominal output file, but the test one sets only
# `trigger_mode = 3`, while run_crossx_hist_filling_pp24.sh also sets `output_generic_hists` and
# `output_gapcut_hists`. Running the test variant here therefore REPLACED the nominal pp24
# histogram file with one missing the entire `_wgapcut` family -- the fiducial-gap histograms this
# analysis is built on -- with no error and no warning.
run_and_log "RDF crossx pp24" "cd '${RDF_DIR}' && bash run_crossx_hist_filling_pp24.sh"

# 2) Plot.  The NOMINAL pair-pT view is ParamsSet::pT_bins_150 (16 log bins 9 -> 150 GeV) and goes
#    to the unsuffixed directories (pp24/, pbpb_..._combined/).  also_pt_120=true additionally
#    refreshes the OPT-IN 9 -> 120 GeV alternative in the *_pt_120 dirs — it must be refreshed in
#    the SAME run as the nominal one, else it silently goes stale, as it did for months.
run_and_log "Plot crossx pp24" "cd '${ANALYSIS_DIR}' && root -l -b -q 'plotting_codes/single_b_analysis/plot_single_b_crossx_pp.cxx(24,\"\",true)'"
run_and_log "Plot crossx pbpb combined" "cd '${ANALYSIS_DIR}' && root -l -b -q 'plotting_codes/single_b_analysis/plot_single_b_crossx_pbpb.cxx(true)'"

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
  if [[ ! "${png}" -nt "${RUN_MARKER}" ]]; then
    echo "[ERR] STALE PNG (not written by this run): ${png}" >&2
    echo "      The plotting macro threw and ROOT exited 0 anyway; this is the previous round's file." >&2
    exit 1
  fi
  echo "[OK ] ${png} ($(stat -c%s "${png}") bytes, fresh)"
}

validate_png "${OUT_DIR}/pp24/pp24_crossx_pair_pt_pair_eta.png"
validate_png "${OUT_DIR}/pp24/pp24_crossx_pair_pt_minv.png"
validate_png "${OUT_DIR}/pp24/pp24_crossx_pair_pt_dr.png"
validate_png "${OUT_DIR}/pp24/pp24_crossx_pair_pt_in_eta_subplots_dr_lines.png"

# Combined PbPb (discovers which years are available; check a representative subset).
# The directory name is built by plot_single_b_crossx_pbpb.cxx::OutDirName() from the years
# it actually found, so derive it here from the SAME list instead of hard-coding a year set
# -- otherwise this validation silently checks a stale directory once a year is added.
COMB_YEARS="$(IFS=_; echo "${PBPB_YEARS[*]}")"
COMB_TAA="${OUT_DIR}/pbpb_${COMB_YEARS}_combined/TAA_weighted"
echo "[INFO] validating combined output in ${COMB_TAA}"
for ctr in ctr0_5 ctr5_10 ctr10_20 ctr20_30 ctr30_50 ctr50_80; do
  validate_png "${COMB_TAA}/pbpb_combined_${ctr}_pair_pt_pair_eta.png"
  validate_png "${COMB_TAA}/pbpb_combined_${ctr}_pair_pt_in_eta_subplots_dr_lines.png"
done

echo "[DONE] All crossx filling and plotting finished successfully."
