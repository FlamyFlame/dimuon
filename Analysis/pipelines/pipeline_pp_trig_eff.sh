#!/usr/bin/env bash
set -Eeuo pipefail

# End-to-end PP trigger efficiency pipeline:
#
# Pipeline 2 (no-correlation single-muon efficiency):
#   1) submit NTuple processing condor jobs (trig eff mode, trigger_mode=1)
#   2) wait for cluster to finish
#   3) validate per-batch ROOT outputs
#   4) hadd per-batch outputs into combined muon_pairs + hists_cut_acceptance
#   5) RDF hist filling — fine q·η (for pT turn-on fitting, TH2D fallback, and plotting)
#   6) pT turn-on fitting (erf+log)
#   7) validate fitting output (TF1s) and TH2D fallback histograms
#   8) Pipeline 2 plotting (no-correlation single-muon trigger efficiency)
#
# PP has NO event selection (no ZDC, no centrality, no nTrk_HITight).
# PP has NO Pipeline 3 (dR corrections excluded — biased procedure).
#
# Usage:
#   ./pipeline_pp_trig_eff.sh
#
# Optional env vars:
#   POLL_SECONDS=45
#   CONDOR_TIMEOUT_SECONDS=0   # 0 => no timeout
#   SKIP_CONDOR=1              # skip condor submit+wait, reuse existing NTuple outputs
#   RDF_NTHREADS=2             # ROOT implicit MT thread count (default 2)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
NTP_DIR="${ANALYSIS_DIR}/NTupleProcessingCode"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
PLOT_DIR="${ANALYSIS_DIR}/plotting_codes/trig_effcy"
FITTER_DIR="${ANALYSIS_DIR}"

POLL_SECONDS="${POLL_SECONDS:-45}"
CONDOR_TIMEOUT_SECONDS="${CONDOR_TIMEOUT_SECONDS:-0}"
SKIP_CONDOR="${SKIP_CONDOR:-0}"
RDF_NTHREADS="${RDF_NTHREADS:-2}"
DATA_BASE="/usatlas/u/yuhanguo/usatlasdata/dimuon_data"
PP_DIR="${DATA_BASE}/pp_2024"

QUEUE_COUNT=12

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

# NTuple output filename pattern for PP trig eff (trigger_mode=1, single_mu4):
#   muon_pairs_pp_2024_partN_single_mu4_mindR_0_02_res_cut_v2.root
#   hists_cut_acceptance_pp_2024_partN_single_mu4_mindR_0_02_res_cut_v2.root
# resonance_cut_mode is forced to 2 when trigger_effcy_calc=true.
TRIG_SUFFIX="_single_mu4_mindR_0_02_res_cut_v2"

get_batch_muon_pairs() {
  local part="$1"
  echo "${PP_DIR}/muon_pairs_pp_2024_part${part}${TRIG_SUFFIX}.root"
}
get_batch_hists() {
  local part="$1"
  echo "${PP_DIR}/hists_cut_acceptance_pp_2024_part${part}${TRIG_SUFFIX}.root"
}

get_combined_muon_pairs() {
  echo "${PP_DIR}/muon_pairs_pp_2024${TRIG_SUFFIX}.root"
}
get_combined_hists() {
  echo "${PP_DIR}/hists_cut_acceptance_pp_2024${TRIG_SUFFIX}.root"
}

get_rdf_output() {
  echo "${PP_DIR}/histograms_real_pairs_pp_2024_single_mu4_fine_q_eta_bin.root"
}
get_fit_output() {
  echo "${PP_DIR}/trg_effcy_pT_fitting_to_erf_plus_log/single_mu_effcy_pT_fit.root"
}

# ==============================
# Pipeline execution
# ==============================
log "PP trigger efficiency pipeline"
log "Sourcing ~/setup.sh"
source_env_once
require_cmd root
require_cmd hadd

if [[ "$SKIP_CONDOR" -eq 1 ]]; then
  log "SKIP_CONDOR=1 — skipping condor submit+wait, reusing existing NTuple outputs"
else
  require_cmd condor_submit
  require_cmd condor_q

  # ------ Stage 1: Submit NTuple trig eff condor jobs ------
  mkdir -p "${NTP_DIR}/logs"
  log "Submitting NTuple trig eff jobs for PP 2024 (trigger_mode=1, single_mu4)"
  pushd "$NTP_DIR" >/dev/null
  out="$(condor_submit run_pp_24.sub)"
  popd >/dev/null
  echo "$out" >&2
  cid="$(extract_cluster_id "$out")"
  log "PP 2024: cluster ${cid}"

  # ------ Stage 2: Wait for cluster ------
  log "Waiting for cluster ${cid}"
  wait_for_cluster_completion "${cid}"
fi

# ------ Stage 3: Validate per-batch outputs ------
declare -a batch_files=()
for (( p=1; p<=QUEUE_COUNT; p++ )); do
  batch_files+=( "$(get_batch_muon_pairs "$p")" )
  batch_files+=( "$(get_batch_hists "$p")" )
done
validate_files_or_fail "NTuple pp24" "${batch_files[@]}"
unset batch_files

# ------ Stage 4: hadd ------
combined_mp="$(get_combined_muon_pairs)"
combined_h="$(get_combined_hists)"

log "hadd muon_pairs for PP 2024"
rm -f "$combined_mp"
declare -a mp_parts=()
for (( p=1; p<=QUEUE_COUNT; p++ )); do mp_parts+=( "$(get_batch_muon_pairs "$p")" ); done
hadd -f "$combined_mp" "${mp_parts[@]}"
validate_files_or_fail "hadd muon_pairs pp24" "$combined_mp"
validate_combined_muon_pair_trees_nonempty_or_fail "$combined_mp"
unset mp_parts

log "hadd hists_cut_acceptance for PP 2024"
rm -f "$combined_h"
declare -a h_parts=()
for (( p=1; p<=QUEUE_COUNT; p++ )); do h_parts+=( "$(get_batch_hists "$p")" ); done
hadd -f "$combined_h" "${h_parts[@]}"
validate_files_or_fail "hadd hists pp24" "$combined_h"
unset h_parts

# ======================================================================
# Pipeline 2: No-correlation single-muon efficiency
# ======================================================================

# ------ Stage 5: RDF Pipeline 2 hist filling (fine q·η) ------
log "Running RDF Pipeline 2 hist filling for PP 2024"
pushd "$RDF_DIR" >/dev/null
root_rc=0
root -l -b <<EOF || root_rc=$?
ROOT::EnableImplicitMT(${RDF_NTHREADS});
.L RDFBasedHistFillingPP.cxx+

RDFBasedHistFillingPP pp(24);
pp.trigger_mode = 1;

std::cout << "\\n[RUN] Starting PP 2024 Pipeline 2 hist filling..." << std::endl;
pp.Run();
std::cout << "\\n[RUN] PP 2024 Pipeline 2 hist filling completed." << std::endl;

gSystem->Exit(0);
EOF
popd >/dev/null
if [[ $root_rc -ne 0 ]]; then
  fail "RDF Pipeline 2 hist filling failed for PP 2024 (ROOT exit code ${root_rc})"
fi
validate_files_or_fail "RDF Pipeline 2 pp24" "$(get_rdf_output)"

# ------ Stage 6: pT turn-on fitting (erf+log) ------
log "Running pT turn-on fitting for PP 2024 (erf+log)"
pushd "$FITTER_DIR" >/dev/null
root_rc=0
root -l -b <<EOF || root_rc=$?
.L SingleMuEffcyPtTurnOnFitter.cxx
single_muon_trig_effcy_pT_fitting();
gSystem->Exit(0);
EOF
popd >/dev/null
if [[ $root_rc -ne 0 ]]; then
  fail "pT turn-on fitting failed for PP 2024 (ROOT exit code ${root_rc})"
fi
validate_files_or_fail "pT fit pp24" "$(get_fit_output)"

# ------ Stage 7: Validate fitting TF1s and TH2D fallback histograms ------
fit_file="$(get_fit_output)"
rdf_file="$(get_rdf_output)"
log "Validating TF1 fits and TH2D fallback for PP 2024"
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
std::cout << "PASS: pp24 TF1=" << ntf1 << " TH2D_fallback=" << n2d_div << std::endl;
gSystem->Exit(0);
EOF
)" || root_rc=$?
echo "$check_output"
if [[ $root_rc -ne 0 ]] || echo "$check_output" | grep -q "^FAIL:"; then
  fail "Fit/TH2D validation failed for PP 2024 (ROOT exit code ${root_rc})"
fi

# ------ Stage 8: Pipeline 2 plotting (no-corr single-muon trig eff) ------
log "Running Pipeline 2 plotting for PP 2024"
pushd "$PLOT_DIR" >/dev/null
root -l -b -q 'trig_effcy_plot_pp.cxx()'
popd >/dev/null

# No Pipeline 3 — dR corrections excluded (biased procedure, same as PbPb).

log "PP trigger efficiency pipeline completed successfully"
