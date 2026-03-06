#!/usr/bin/env bash
set -Eeuo pipefail

# End-to-end pipeline:
#   1) submit CRTP condor jobs for selected mode
#   2) wait for NEWLY submitted cluster to finish
#   3) validate per-batch ROOT outputs
#   4) hadd -> combined muon_pairs input for RDF
#   5) run RDF
#   6) validate RDF outputs (nominal + near/away divided)
#   7) run flavor/origin categorized plotting
#
# Usage:
#   ./pipeline_pythia_truth.sh private
#   ./pipeline_pythia_truth.sh nonprivate_5p36
#   ./pipeline_pythia_truth.sh nonprivate_5p02
#
# Optional env vars:
#   POLL_SECONDS=45
#   CONDOR_TIMEOUT_SECONDS=0   # 0 => no timeout

MODE="${1:-}"
if [[ -z "${MODE}" ]]; then
  echo "Usage: $0 [private|nonprivate_5p36|nonprivate_5p02]"
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
PLOT_DIR="${ANALYSIS_DIR}/plotting_codes/pythia_plotting_codes"

POLL_SECONDS="${POLL_SECONDS:-45}"
CONDOR_TIMEOUT_SECONDS="${CONDOR_TIMEOUT_SECONDS:-0}"

now() { date '+%F %T'; }
log() { echo "[$(now)] $*"; }

on_error() {
  local exit_code=$?
  echo "[$(now)] ERROR: command failed (exit=${exit_code}) at line ${BASH_LINENO[0]}: ${BASH_COMMAND}" >&2
}
trap on_error ERR

fail() {
  echo "[$(now)] ERROR: $*" >&2
  exit 1
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || fail "Missing required command: $1"
}

source_env_once() {
  local err_trap_saved
  err_trap_saved="$(trap -p ERR || true)"
  trap - ERR

  set +e
  set +u
  source ~/setup.sh
  local rc=$?
  set -e
  set -u

  if [[ -n "${err_trap_saved}" ]]; then
    eval "${err_trap_saved}"
  else
    trap on_error ERR
  fi

  if [[ $rc -ne 0 ]]; then
    fail "Failed to source ~/setup.sh"
  fi
}

validate_root_file_quick() {
  local f="$1"
  [[ -f "$f" ]] || return 1
  [[ -s "$f" ]] || return 1

  root -l -b -q <<EOF >/dev/null 2>&1
TFile *fin = TFile::Open("$f", "READ");
if (!fin || fin->IsZombie()) { gSystem->Exit(2); }
if (!fin->GetListOfKeys() || fin->GetListOfKeys()->GetSize() <= 0) { fin->Close(); gSystem->Exit(3); }
fin->Close();
gSystem->Exit(0);
EOF
  return $?
}

validate_files_or_fail() {
  local label="$1"
  shift
  local files=("$@")
  local bad=0
  for f in "${files[@]}"; do
    if validate_root_file_quick "$f"; then
      log "OK ${label}: $f"
    else
      echo "[$(now)] BAD ${label}: $f" >&2
      bad=1
    fi
  done
  [[ $bad -eq 0 ]] || fail "Validation failed for ${label} files"
}

validate_combined_muon_pair_trees_nonempty_or_fail() {
  local f="$1"
  [[ -f "$f" ]] || fail "Combined file not found: $f"

  local check_output
  if check_output="$(root -l -b -q <<EOF
TFile *fin = TFile::Open("$f", "READ");
if (!fin || fin->IsZombie()) {
  std::cout << "ERROR: cannot open file" << std::endl;
  gSystem->Exit(2);
}
TTree *t1 = dynamic_cast<TTree*>(fin->Get("muon_pair_tree_sign1"));
TTree *t2 = dynamic_cast<TTree*>(fin->Get("muon_pair_tree_sign2"));
if (!t1 || !t2) {
  std::cout << "ERROR: missing required trees (muon_pair_tree_sign1 and/or muon_pair_tree_sign2)" << std::endl;
  fin->Close();
  gSystem->Exit(3);
}
Long64_t n1 = t1->GetEntries();
Long64_t n2 = t2->GetEntries();
std::cout << "Combined tree entries: muon_pair_tree_sign1=" << n1
          << ", muon_pair_tree_sign2=" << n2 << std::endl;
if (n1 <= 0 || n2 <= 0) {
  std::cout << "ERROR: combined muon_pairs file has empty sign tree(s)" << std::endl;
  fin->Close();
  gSystem->Exit(4);
}
fin->Close();
gSystem->Exit(0);
EOF
)"; then
    echo "$check_output"
  else
    echo "$check_output"
    fail "Post-hadd tree non-empty validation failed for: $f"
  fi
}

extract_cluster_id() {
  local submit_output="$1"
  local cid
  cid="$(echo "$submit_output" | sed -n 's/.*cluster \([0-9]\+\).*/\1/p' | tail -n1)"
  [[ -n "$cid" ]] || fail "Could not parse ClusterId from condor_submit output"
  echo "$cid"
}

wait_for_cluster_completion() {
  local cluster_id="$1"
  local t0
  t0="$(date +%s)"

  while true; do
    local qout
    qout="$(condor_q "${cluster_id}" -autoformat ClusterId ProcId JobStatus 2>/dev/null || true)"

    if [[ -z "$qout" ]]; then
      log "Cluster ${cluster_id} no longer in condor_q (finished/left queue)."
      break
    fi

    local total=0 idle=0 running=0 held=0 other=0
    while read -r cid pid st; do
      [[ -n "${cid:-}" ]] || continue
      total=$((total + 1))
      case "$st" in
        1) idle=$((idle + 1)) ;;
        2) running=$((running + 1)) ;;
        5) held=$((held + 1)) ;;
        *) other=$((other + 1)) ;;
      esac
    done <<< "$qout"

    log "Cluster ${cluster_id}: queued=${total}, idle=${idle}, running=${running}, held=${held}, other=${other}"

    if (( held > 0 )); then
      fail "Cluster ${cluster_id} has held jobs; inspect with: condor_q ${cluster_id} -hold"
    fi

    if (( CONDOR_TIMEOUT_SECONDS > 0 )); then
      local tnow elapsed
      tnow="$(date +%s)"
      elapsed=$((tnow - t0))
      if (( elapsed > CONDOR_TIMEOUT_SECONDS )); then
        fail "Timed out waiting for cluster ${cluster_id} after ${elapsed} s"
      fi
    fi

    sleep "$POLL_SECONDS"
  done
}

run_plotting_stage() {
  local is_private_root="$1"
  local ecom="$2"

  local flavor_macro
  local origin_macro
  flavor_macro="${PLOT_DIR}/plot_flavor_categorized_kinematics.cxx(${is_private_root},${ecom})"
  origin_macro="${PLOT_DIR}/plot_origin_categorized_kinematics.cxx(${is_private_root},${ecom})"

  log "Running flavor categorized plotting"
  root -l -b -q "$flavor_macro"

  log "Running origin categorized plotting"
  root -l -b -q "$origin_macro"
}

run_rdf_stage() {
  local mode="$1"      # private|nonprivate
  local ecom="$2"      # 5.36|5.02

  pushd "$RDF_DIR" >/dev/null
  log "Running RDF stage: mode=${mode}, ecom=${ecom}, cuts=nocuts"
  ./run_rdf_pythia_truth.sh "$mode" "$ecom" nocuts
  popd >/dev/null
}

submit_condor_cluster() {
  local submit_file="$1"
  local submit_args="$2"
  pushd "$SCRIPT_DIR" >/dev/null
  local out
  if [[ -n "$submit_args" ]]; then
    out="$(condor_submit -append "arguments = ${submit_args}" "$submit_file")"
  else
    out="$(condor_submit "$submit_file")"
  fi
  popd >/dev/null

  echo "$out" >&2
  extract_cluster_id "$out"
}

# ------------------------------
# Mode-specific configuration
# ------------------------------

declare -a BATCH_OUTPUTS
MODE_LABEL=""
CRTP_SUBMIT_FILE=""
CRTP_SUBMIT_ARGS=""
RDF_MODE=""
RDF_ECOM=""
PLOT_IS_PRIVATE=""
PLOT_ECOM=""
COMBINED_OUTPUT=""
HADD_PATTERN=""
RDF_NOMINAL=""
RDF_DIVIDED=""

case "$MODE" in
  private)
    MODE_LABEL="private"
    CRTP_SUBMIT_FILE="run_pythia_truth_kn_private.sub"
    CRTP_SUBMIT_ARGS='$(Process) 1'
    RDF_MODE="private"
    RDF_ECOM="5.36"
    PLOT_IS_PRIVATE="true"
    PLOT_ECOM="5.36"

    BASE_DIR="/usatlas/u/yuhanguo/usatlasdata/pythia_private_sample"
    BATCH_OUTPUTS=(
      "${BASE_DIR}/muon_pairs_pythia_allto0318_no_data_resonance_cuts.root"
      "${BASE_DIR}/muon_pairs_pythia_after0322_no_data_resonance_cuts.root"
    )
    COMBINED_OUTPUT="${BASE_DIR}/muon_pairs_pythia_combined_no_data_resonance_cuts.root"
    HADD_PATTERN="${BASE_DIR}/muon_pairs_pythia_*_no_data_resonance_cuts.root"
    RDF_NOMINAL="${BASE_DIR}/histograms_pythia_combined_no_data_resonance_cuts.root"
    RDF_DIVIDED="${BASE_DIR}/histograms_pythia_combined_no_data_resonance_cuts_near_away_divided.root"
    ;;

  nonprivate_5p36)
    MODE_LABEL="nonprivate_5p36"
    CRTP_SUBMIT_FILE="run_pythia_truth_kn_nonprivate.sub"
    CRTP_SUBMIT_ARGS='$(Process) 0 5.36'
    RDF_MODE="nonprivate"
    RDF_ECOM="5.36"
    PLOT_IS_PRIVATE="false"
    PLOT_ECOM="5.36"

    BASE_DIR="/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5p36TeV"
    BATCH_OUTPUTS=(
      "${BASE_DIR}/muon_pairs_pythia_5p36TeV_kn0_no_data_resonance_cuts.root"
      "${BASE_DIR}/muon_pairs_pythia_5p36TeV_kn1_no_data_resonance_cuts.root"
      "${BASE_DIR}/muon_pairs_pythia_5p36TeV_kn2_no_data_resonance_cuts.root"
      "${BASE_DIR}/muon_pairs_pythia_5p36TeV_kn3_no_data_resonance_cuts.root"
      "${BASE_DIR}/muon_pairs_pythia_5p36TeV_kn4_no_data_resonance_cuts.root"
      "${BASE_DIR}/muon_pairs_pythia_5p36TeV_kn5_no_data_resonance_cuts.root"
    )
    COMBINED_OUTPUT="${BASE_DIR}/muon_pairs_pythia_5p36TeV_no_data_resonance_cuts.root"
    HADD_PATTERN="${BASE_DIR}/muon_pairs_pythia_5p36TeV_kn*_no_data_resonance_cuts.root"
    RDF_NOMINAL="${BASE_DIR}/histograms_pythia_5p36TeV_no_data_resonance_cuts.root"
    RDF_DIVIDED="${BASE_DIR}/histograms_pythia_5p36TeV_no_data_resonance_cuts_near_away_divided.root"
    ;;

  nonprivate_5p02)
    MODE_LABEL="nonprivate_5p02"
    CRTP_SUBMIT_FILE="run_pythia_truth_kn_nonprivate.sub"
    CRTP_SUBMIT_ARGS='$(Process) 0 5.02'
    RDF_MODE="nonprivate"
    RDF_ECOM="5.02"
    PLOT_IS_PRIVATE="false"
    PLOT_ECOM="5.02"

    BASE_DIR="/usatlas/u/yuhanguo/usatlasdata/pythia_truth_full_sample/pythia_5TeV"
    BATCH_OUTPUTS=(
      "${BASE_DIR}/muon_pairs_pythia_5p02TeV_kn0_no_data_resonance_cuts.root"
      "${BASE_DIR}/muon_pairs_pythia_5p02TeV_kn1_no_data_resonance_cuts.root"
      "${BASE_DIR}/muon_pairs_pythia_5p02TeV_kn2_no_data_resonance_cuts.root"
      "${BASE_DIR}/muon_pairs_pythia_5p02TeV_kn3_no_data_resonance_cuts.root"
      "${BASE_DIR}/muon_pairs_pythia_5p02TeV_kn4_no_data_resonance_cuts.root"
      "${BASE_DIR}/muon_pairs_pythia_5p02TeV_kn5_no_data_resonance_cuts.root"
    )
    COMBINED_OUTPUT="${BASE_DIR}/muon_pairs_pythia_5p02TeV_no_data_resonance_cuts.root"
    HADD_PATTERN="${BASE_DIR}/muon_pairs_pythia_5p02TeV_kn*_no_data_resonance_cuts.root"
    RDF_NOMINAL="${BASE_DIR}/histograms_pythia_5p02TeV_no_data_resonance_cuts.root"
    RDF_DIVIDED="${BASE_DIR}/histograms_pythia_5p02TeV_no_data_resonance_cuts_near_away_divided.root"
    ;;

  *)
    fail "Unknown mode '${MODE}'. Use private|nonprivate_5p36|nonprivate_5p02"
    ;;
esac

# ------------------------------
# Pipeline execution
# ------------------------------
require_cmd condor_submit
require_cmd condor_q

log "Mode: ${MODE_LABEL}"
log "Sourcing ~/setup.sh for ROOT/LCG tools (including hadd)"
source_env_once
require_cmd root
require_cmd hadd

log "Submitting CRTP condor jobs via ${CRTP_SUBMIT_FILE}"
cluster_id="$(submit_condor_cluster "${CRTP_SUBMIT_FILE}" "${CRTP_SUBMIT_ARGS}")"
log "Submitted cluster: ${cluster_id}"

log "Waiting for completion of submitted cluster only"
wait_for_cluster_completion "${cluster_id}"

log "Validating per-batch CRTP outputs before hadd"
validate_files_or_fail "CRTP batch" "${BATCH_OUTPUTS[@]}"

log "Creating combined CRTP input for RDF with hadd"
rm -f "${COMBINED_OUTPUT}"
# shellcheck disable=SC2086
hadd -f "${COMBINED_OUTPUT}" ${HADD_PATTERN}
validate_files_or_fail "CRTP combined" "${COMBINED_OUTPUT}"
validate_combined_muon_pair_trees_nonempty_or_fail "${COMBINED_OUTPUT}"

run_rdf_stage "${RDF_MODE}" "${RDF_ECOM}"

log "Validating RDF outputs"
validate_files_or_fail "RDF" "${RDF_NOMINAL}" "${RDF_DIVIDED}"

run_plotting_stage "${PLOT_IS_PRIVATE}" "${PLOT_ECOM}"

log "Pipeline completed successfully for mode=${MODE_LABEL}"
