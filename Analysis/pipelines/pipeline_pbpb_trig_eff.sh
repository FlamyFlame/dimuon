#!/usr/bin/env bash
set -Eeuo pipefail

# End-to-end PbPb trigger efficiency pipeline:
#
# Event selection (shared prerequisite):
#   0) derive event selection cuts + produce all event selection plots
#
# Pipeline 2 (no-correlation single-muon efficiency):
#   1) submit NTuple processing condor jobs (trig eff mode, all 3 PbPb years)
#   2) wait for all clusters to finish
#   3) validate per-batch ROOT outputs
#   4) hadd per-batch outputs into combined muon_pairs + hists_cut_acceptance per year
#   5) RDF hist filling — fine q·η (for pT turn-on fitting, TH2D fallback, and plotting)
#   6) pT turn-on fitting (Fermi+log) per year
#   7) validate fitting output (TF1s) and TH2D fallback histograms
#   8) Pipeline 2 plotting (no-correlation single-muon trigger efficiency)
#
# Pipeline 3 (dR corrections via inverse weighting):
#   9) RDF hist filling — inv_weight_by_single_mu_effcy cycle per year
#   10) Pipeline 3 plotting (dR correction overlays + plateau normalization)
#
# Usage:
#   ./pipeline_pbpb_trig_eff.sh
#
# Optional env vars:
#   POLL_SECONDS=45
#   CONDOR_TIMEOUT_SECONDS=0   # 0 => no timeout
#   YEARS="23 24 25"           # override which years to process
#   SKIP_CONDOR=1              # skip condor submit+wait, reuse existing NTuple outputs
#   SKIP_EVSEL=1               # skip event selection only, still run condor+RDF+plots
#   RDF_NTHREADS=2             # ROOT implicit MT thread count (default 2; 8 needs ~48 GB RAM)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
NTP_DIR="${ANALYSIS_DIR}/NTupleProcessingCode"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
PLOT_DIR="${ANALYSIS_DIR}/plotting_codes/trig_effcy"
EVSEL_DIR="${ANALYSIS_DIR}/plotting_codes/event_selection"

POLL_SECONDS="${POLL_SECONDS:-45}"
CONDOR_TIMEOUT_SECONDS="${CONDOR_TIMEOUT_SECONDS:-0}"
SKIP_CONDOR="${SKIP_CONDOR:-0}"
SKIP_EVSEL="${SKIP_EVSEL:-${SKIP_CONDOR}}"
RDF_NTHREADS="${RDF_NTHREADS:-2}"
YEARS=(${YEARS:-23 24 25})
DATA_BASE="/usatlas/u/yuhanguo/usatlasdata/dimuon_data"
PLOT_BASE="${DATA_BASE}/plots"
FITTER_DIR="${ANALYSIS_DIR}"
PLOT_DR_DIR="${ANALYSIS_DIR}/plotting_codes"

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
  set +e; set +u
  source ~/setup.sh
  local rc=$?
  set -e; set -u
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
  local label="$1"; shift
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
  std::cout << "ERROR: missing required trees" << std::endl;
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
  local t0; t0="$(date +%s)"
  while true; do
    local qout
    qout="$(condor_q "${cluster_id}" -autoformat ClusterId ProcId JobStatus 2>/dev/null || true)"
    if [[ -z "$qout" ]]; then
      log "Cluster ${cluster_id} finished."
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
      tnow="$(date +%s)"; elapsed=$((tnow - t0))
      if (( elapsed > CONDOR_TIMEOUT_SECONDS )); then
        fail "Timed out waiting for cluster ${cluster_id} after ${elapsed} s"
      fi
    fi
    sleep "$POLL_SECONDS"
  done
}

# Queue counts per year (must match .sub files)
declare -A QUEUE_COUNTS=( [23]=4 [24]=2 [25]=6 )

get_year_dir() { echo "${DATA_BASE}/pbpb_20$1"; }

# NTuple output filename pattern for trig eff mode (default):
#   muon_pairs_pbpb_20YY_partN_single_mu4_mindR_0_02_res_cut_v2.root
#   hists_cut_acceptance_pbpb_20YY_partN_single_mu4_mindR_0_02_res_cut_v2.root
get_batch_muon_pairs() {
  local yr="$1" part="$2"
  echo "$(get_year_dir "$yr")/muon_pairs_pbpb_20${yr}_part${part}_single_mu4_mindR_0_02_res_cut_v2.root"
}
get_batch_hists() {
  local yr="$1" part="$2"
  echo "$(get_year_dir "$yr")/hists_cut_acceptance_pbpb_20${yr}_part${part}_single_mu4_mindR_0_02_res_cut_v2.root"
}

# Combined (hadd) output filenames
get_combined_muon_pairs() {
  echo "$(get_year_dir "$1")/muon_pairs_pbpb_20${1}_single_mu4_mindR_0_02_res_cut_v2.root"
}
get_combined_hists() {
  echo "$(get_year_dir "$1")/hists_cut_acceptance_pbpb_20${1}_single_mu4_mindR_0_02_res_cut_v2.root"
}

# RDF output filename (fine q·η binning — needed for pT fitting + TH2D fallback)
get_rdf_output() {
  echo "$(get_year_dir "$1")/histograms_real_pairs_pbpb_20${1}_single_mu4_fine_q_eta_bin.root"
}
get_fit_output() {
  echo "$(get_year_dir "$1")/trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit.root"
}

# ==============================
# Pipeline execution
# ==============================
log "PbPb trigger efficiency pipeline — years: ${YEARS[*]}"
log "Sourcing ~/setup.sh"
source_env_once
require_cmd root
require_cmd hadd

# ------ Stage 0: Event selection (derive cuts + plots) ------
if [[ "$SKIP_EVSEL" -eq 1 ]]; then
  log "SKIP_EVSEL=1 — skipping event selection (reusing existing cuts)"
else
  log "Running event selection for all years"
  for yr in "${YEARS[@]}"; do
    log "Event selection: deriving cuts for PbPb 20${yr}"
    pushd "$EVSEL_DIR" >/dev/null
    root -l -b -q "plot_pbpb_event_sel_event_level.cxx(${yr})"
    popd >/dev/null
    for suffix in "" "_alt"; do
      evsel_file="${DATA_BASE}/pbpb_20${yr}/event_sel_cuts_pbpb_20${yr}${suffix}.root"
      validate_root_file_quick "$evsel_file" || fail "Event selection output missing: $evsel_file"
      log "OK event_sel yr${yr}: $evsel_file"
    done

    log "Event selection: plotting cuts for PbPb 20${yr} (nominal + alt)"
    pushd "$EVSEL_DIR" >/dev/null
    root -l -b -q "plot_pbpb_event_sel_cuts.cxx(${yr})"
    root -l -b <<EOF
.L plot_pbpb_event_sel_cuts.cxx
plot_pbpb_event_sel_cuts_alt(${yr});
gSystem->Exit(0);
EOF
    popd >/dev/null
  done

  log "Event selection: FCal comparison plots (all years)"
  pushd "$EVSEL_DIR" >/dev/null
  root -l -b -q 'plot_pbpb_fcal_comparison.cxx()'
  root -l -b <<EOF
.L plot_pbpb_fcal_comparison.cxx
plot_pbpb_fcal_comparison_alt();
gSystem->Exit(0);
EOF
  root -l -b <<EOF
.L plot_pbpb_fcal_comparison.cxx
plot_pbpb_fcal_comparison_2524();
gSystem->Exit(0);
EOF
  root -l -b <<EOF
.L plot_pbpb_fcal_comparison.cxx
plot_pbpb_fcal_comparison_2524_alt();
gSystem->Exit(0);
EOF
  popd >/dev/null
  log "Event selection complete for all years"
fi

if [[ "$SKIP_CONDOR" -eq 1 ]]; then
  log "SKIP_CONDOR=1 — skipping condor submit+wait, reusing existing NTuple outputs"
else
  require_cmd condor_submit
  require_cmd condor_q

  # ------ Stage 1: Submit NTuple condor jobs ------
  mkdir -p "${NTP_DIR}/logs"
  declare -A CLUSTER_IDS
  for yr in "${YEARS[@]}"; do
    log "Submitting NTuple trig eff jobs for PbPb 20${yr}"
    pushd "$NTP_DIR" >/dev/null
    out="$(condor_submit "run_pbpb_${yr}.sub")"
    popd >/dev/null
    echo "$out" >&2
    cid="$(extract_cluster_id "$out")"
    CLUSTER_IDS[$yr]="$cid"
    log "Year ${yr}: cluster ${cid}"
  done

  # ------ Stage 2: Wait for all clusters ------
  for yr in "${YEARS[@]}"; do
    log "Waiting for year ${yr} cluster ${CLUSTER_IDS[$yr]}"
    wait_for_cluster_completion "${CLUSTER_IDS[$yr]}"
  done
fi

# ------ Stage 3: Validate per-batch outputs ------
for yr in "${YEARS[@]}"; do
  nq="${QUEUE_COUNTS[$yr]}"
  declare -a batch_files=()
  for (( p=1; p<=nq; p++ )); do
    batch_files+=( "$(get_batch_muon_pairs "$yr" "$p")" )
    batch_files+=( "$(get_batch_hists "$yr" "$p")" )
  done
  validate_files_or_fail "NTuple yr${yr}" "${batch_files[@]}"
  unset batch_files
done

# ------ Stage 4: hadd per year ------
for yr in "${YEARS[@]}"; do
  nq="${QUEUE_COUNTS[$yr]}"
  combined_mp="$(get_combined_muon_pairs "$yr")"
  combined_h="$(get_combined_hists "$yr")"

  log "hadd muon_pairs for year ${yr}"
  rm -f "$combined_mp"
  declare -a mp_parts=()
  for (( p=1; p<=nq; p++ )); do mp_parts+=( "$(get_batch_muon_pairs "$yr" "$p")" ); done
  hadd -f "$combined_mp" "${mp_parts[@]}"
  validate_files_or_fail "hadd muon_pairs yr${yr}" "$combined_mp"
  validate_combined_muon_pair_trees_nonempty_or_fail "$combined_mp"
  unset mp_parts

  log "hadd hists_cut_acceptance for year ${yr}"
  rm -f "$combined_h"
  declare -a h_parts=()
  for (( p=1; p<=nq; p++ )); do h_parts+=( "$(get_batch_hists "$yr" "$p")" ); done
  hadd -f "$combined_h" "${h_parts[@]}"
  validate_files_or_fail "hadd hists yr${yr}" "$combined_h"
  unset h_parts
done

# ======================================================================
# Pipeline 2: No-correlation single-muon efficiency
# ======================================================================

# ------ Stage 5: RDF Pipeline 2 hist filling (fine q·η) ------
for yr in "${YEARS[@]}"; do
  log "Running RDF Pipeline 2 hist filling for PbPb 20${yr}"
  pushd "$RDF_DIR" >/dev/null
  root_rc=0
  root -l -b <<EOF || root_rc=$?
ROOT::EnableImplicitMT(${RDF_NTHREADS});
.L RDFBasedHistFillingPbPb.cxx+

RDFBasedHistFillingPbPb pbpb(${yr});
pbpb.trigger_mode = 1;
pbpb.mindR_trig = 0.02;

std::cout << "\\n[RUN] Starting PbPb 20${yr} Pipeline 2 hist filling..." << std::endl;
pbpb.Run();
std::cout << "\\n[RUN] PbPb 20${yr} Pipeline 2 hist filling completed." << std::endl;

gSystem->Exit(0);
EOF
  popd >/dev/null
  if [[ $root_rc -ne 0 ]]; then
    fail "RDF Pipeline 2 hist filling failed for year ${yr} (ROOT exit code ${root_rc})"
  fi
  validate_files_or_fail "RDF Pipeline 2 yr${yr}" "$(get_rdf_output "$yr")"
done

# ------ Stage 6: pT turn-on fitting (Fermi+log) per year ------
# Must run interpreted (.L without +) due to CRTP rootcling issues
for yr in "${YEARS[@]}"; do
  log "Running pT turn-on fitting for PbPb 20${yr}"
  pushd "$FITTER_DIR" >/dev/null
  root_rc=0
  root -l -b <<EOF || root_rc=$?
.L SingleMuEffcyPtTurnOnFitter.cxx
single_muon_trig_effcy_pT_fitting_PbPb(${yr});
gSystem->Exit(0);
EOF
  popd >/dev/null
  if [[ $root_rc -ne 0 ]]; then
    fail "pT turn-on fitting failed for year ${yr} (ROOT exit code ${root_rc})"
  fi
  validate_files_or_fail "pT fit yr${yr}" "$(get_fit_output "$yr")"
done

# ------ Stage 7: Validate fitting TF1s and TH2D fallback histograms ------
for yr in "${YEARS[@]}"; do
  fit_file="$(get_fit_output "$yr")"
  rdf_file="$(get_rdf_output "$yr")"
  log "Validating TF1 fits and TH2D fallback for PbPb 20${yr}"
  root_rc=0
  check_output="$(root -l -b -q <<EOF
TFile *ffit = TFile::Open("${fit_file}", "READ");
if (!ffit || ffit->IsZombie()) {
  std::cout << "FAIL: cannot open fit file ${fit_file}" << std::endl;
  gSystem->Exit(1);
}
int ntf1 = 0;
TIter next_fit(ffit->GetListOfKeys());
TKey *k;
while ((k = (TKey*)next_fit())) {
  if (std::string(k->GetClassName()) == "TF1") ntf1++;
}
std::cout << "TF1 count: " << ntf1 << std::endl;
if (ntf1 < 10) {
  std::cout << "FAIL: expected at least 10 TF1s, found " << ntf1 << std::endl;
  ffit->Close();
  gSystem->Exit(2);
}
ffit->Close();

TFile *fhist = TFile::Open("${rdf_file}", "READ");
if (!fhist || fhist->IsZombie()) {
  std::cout << "FAIL: cannot open hist file ${rdf_file}" << std::endl;
  gSystem->Exit(3);
}
int n2d_div = 0;
TIter next_h(fhist->GetListOfKeys());
while ((k = (TKey*)next_h())) {
  std::string name = k->GetName();
  if (std::string(k->GetClassName()) == "TH2D"
      && name.find("_divided") != std::string::npos
      && name.find("pt2nd_vs_q_eta2nd") != std::string::npos) {
    n2d_div++;
  }
}
std::cout << "TH2D _divided (pt2nd_vs_q_eta2nd) count: " << n2d_div << std::endl;
if (n2d_div < 4) {
  std::cout << "FAIL: expected at least 4 TH2D _divided histograms, found " << n2d_div << std::endl;
  fhist->Close();
  gSystem->Exit(4);
}
fhist->Close();
std::cout << "PASS: yr${yr} TF1=" << ntf1 << " TH2D_fallback=" << n2d_div << std::endl;
gSystem->Exit(0);
EOF
)" || root_rc=$?
  echo "$check_output"
  if [[ $root_rc -ne 0 ]] || echo "$check_output" | grep -q "^FAIL:"; then
    fail "Fit/TH2D validation failed for year ${yr} (ROOT exit code ${root_rc})"
  fi
done

# ------ Stage 8: Pipeline 2 plotting (no-corr single-muon trig eff) per year ------
for yr in "${YEARS[@]}"; do
  log "Running Pipeline 2 plotting for PbPb 20${yr}"
  pushd "$PLOT_DIR" >/dev/null
  root -l -b -q "trig_effcy_plot_PbPb.cxx(${yr})"
  popd >/dev/null
done

# ======================================================================
# Pipeline 3: dR corrections via inverse weighting
# ======================================================================

# ------ Stage 9: RDF Pipeline 3 hist filling (inv_weight_by_single_mu_effcy) ------
for yr in "${YEARS[@]}"; do
  log "Running RDF Pipeline 3 hist filling for PbPb 20${yr}"
  pushd "$RDF_DIR" >/dev/null
  root_rc=0
  root -l -b <<EOF || root_rc=$?
ROOT::EnableImplicitMT(${RDF_NTHREADS});
.L RDFBasedHistFillingPbPb.cxx+

RDFBasedHistFillingPbPb pbpb(${yr});
pbpb.trigger_mode = 1;
pbpb.mindR_trig = 0.02;
pbpb.hist_filling_cycle = 2;  // inv_weight_by_single_mu_effcy

std::cout << "\\n[RUN] Starting PbPb 20${yr} Pipeline 3 hist filling..." << std::endl;
pbpb.Run();
std::cout << "\\n[RUN] PbPb 20${yr} Pipeline 3 hist filling completed." << std::endl;

gSystem->Exit(0);
EOF
  popd >/dev/null
  if [[ $root_rc -ne 0 ]]; then
    fail "RDF Pipeline 3 hist filling failed for year ${yr} (ROOT exit code ${root_rc})"
  fi
  validate_files_or_fail "RDF Pipeline 3 yr${yr}" "$(get_rdf_output "$yr")"
done

# ------ Stage 10: Pipeline 3 plotting (dR corrections + plateau normalization) ------
log "Running dR correction plots (all years combined)"
pushd "$PLOT_DR_DIR" >/dev/null
root -l -b -q "plot_dR_trig_corr.C"
popd >/dev/null

log "Running plateau normalization (all years combined)"
pushd "$PLOT_DR_DIR" >/dev/null
root -l -b -q "plateau_normalize_dR_corr.C"
popd >/dev/null

log "PbPb trigger efficiency pipeline completed successfully for years: ${YEARS[*]}"
