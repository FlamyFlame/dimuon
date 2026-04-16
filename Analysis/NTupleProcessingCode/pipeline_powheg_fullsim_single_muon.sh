#!/usr/bin/env bash
# End-to-end single-muon reconstruction efficiency pipeline:
#
#   1) Submit Condor jobs to (re)generate single-muon ROOT trees for bb + cc
#   2) Wait for both clusters to finish
#   3) Validate per-batch output files exist and muon_tree is non-empty
#   4) hadd per-batch files into the per-flavour combined trees expected by
#      the mixer and the RDF hist filler
#   5) Run RDF histogram filling for single-muon reco efficiency / det response
#   6) Validate the histogram output file
#   7) Run the plotter for both medium and tight working points
#
# Prerequisites:
#   - PowhegFullSimExtras.c has the truth_pt > 4 && |eta| < 2.4 guard in place
#   - ~/setup.sh provides ROOT, hadd, condor_submit, condor_q
#
# Optional env vars:
#   POLL_SECONDS=45
#   CONDOR_TIMEOUT_SECONDS=0   # 0 => no timeout

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"

POLL_SECONDS="${POLL_SECONDS:-45}"
CONDOR_TIMEOUT_SECONDS="${CONDOR_TIMEOUT_SECONDS:-0}"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
POWHEG_BASE="/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample"
BB_DIR="${POWHEG_BASE}/user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp17.bb.Feb2026.v1._MYSTREAM"
CC_DIR="${POWHEG_BASE}/user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp17.cc.Feb2026.v1._MYSTREAM"

BB_COMBINED="${BB_DIR}/single_muon_trees_powheg_bb_fullsim_pp17.root"
CC_COMBINED="${CC_DIR}/single_muon_trees_powheg_cc_fullsim_pp17.root"

SINGLE_MUON_HIST="${POWHEG_BASE}/histograms_powheg_fullsim_single_muon_pp17.root"

# Build arrays of per-batch output files (file_batch = process + 1, queue 11)
BB_BATCH_FILES=()
CC_BATCH_FILES=()
for batch in $(seq 1 11); do
    BB_BATCH_FILES+=("${BB_DIR}/single_muon_trees_powheg_bb_fullsim_pp17_part${batch}.root")
    CC_BATCH_FILES+=("${CC_DIR}/single_muon_trees_powheg_cc_fullsim_pp17_part${batch}.root")
done

# ---------------------------------------------------------------------------
# Helpers (modelled after pipeline_pythia_truth.sh)
# ---------------------------------------------------------------------------
now()  { date '+%F %T'; }
log()  { echo "[$(now)] $*"; }
fail() { echo "[$(now)] ERROR: $*" >&2; exit 1; }

on_error() {
    local exit_code=$?
    echo "[$(now)] ERROR: command failed (exit=${exit_code}) at line ${BASH_LINENO[0]}: ${BASH_COMMAND}" >&2
}
trap on_error ERR

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
    if [[ -n "${err_trap_saved}" ]]; then eval "${err_trap_saved}"; else trap on_error ERR; fi
    [[ $rc -eq 0 ]] || fail "Failed to source ~/setup.sh"
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
    local bad=0
    for f in "$@"; do
        if validate_root_file_quick "$f"; then
            log "OK ${label}: $f"
        else
            echo "[$(now)] BAD ${label}: $f" >&2
            bad=1
        fi
    done
    [[ $bad -eq 0 ]] || fail "Validation failed for ${label} files"
}

validate_single_muon_tree_nonempty_or_fail() {
    local f="$1"
    [[ -f "$f" ]] || fail "Single-muon combined file not found: $f"

    local check_output
    if check_output="$(root -l -b -q <<EOF
TFile *fin = TFile::Open("$f", "READ");
if (!fin || fin->IsZombie()) {
  std::cout << "ERROR: cannot open file" << std::endl;
  gSystem->Exit(2);
}
TTree *t = dynamic_cast<TTree*>(fin->Get("muon_tree"));
if (!t) {
  std::cout << "ERROR: missing tree muon_tree in $f" << std::endl;
  fin->Close();
  gSystem->Exit(3);
}
Long64_t n = t->GetEntries();
std::cout << "muon_tree entries: " << n << std::endl;
if (n <= 0) {
  std::cout << "ERROR: muon_tree is empty in $f" << std::endl;
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
        fail "muon_tree non-empty validation failed for: $f"
    fi
}

# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------
require_cmd condor_submit
require_cmd condor_q

log "Sourcing ~/setup.sh for ROOT / hadd / condor tools"
source_env_once
require_cmd root
require_cmd hadd

# --- Step 1: submit bb and cc condor clusters ---
log "Submitting bb single-muon Condor jobs (run_powheg_fullsim_wtruth_bb_single_muon.sub)"
pushd "${SCRIPT_DIR}" >/dev/null
bb_submit_out="$(condor_submit run_powheg_fullsim_wtruth_bb_single_muon.sub)"
echo "${bb_submit_out}" >&2
bb_cluster="$(extract_cluster_id "${bb_submit_out}")"
log "Submitted bb cluster: ${bb_cluster}"

log "Submitting cc single-muon Condor jobs (run_powheg_fullsim_wtruth_cc_single_muon.sub)"
cc_submit_out="$(condor_submit run_powheg_fullsim_wtruth_cc_single_muon.sub)"
echo "${cc_submit_out}" >&2
cc_cluster="$(extract_cluster_id "${cc_submit_out}")"
log "Submitted cc cluster: ${cc_cluster}"
popd >/dev/null

# --- Step 2: wait for both clusters ---
log "Waiting for bb cluster ${bb_cluster} to finish"
wait_for_cluster_completion "${bb_cluster}"

log "Waiting for cc cluster ${cc_cluster} to finish"
wait_for_cluster_completion "${cc_cluster}"

# --- Step 3: validate per-batch outputs ---
log "Validating bb per-batch outputs"
validate_files_or_fail "bb batch" "${BB_BATCH_FILES[@]}"

log "Validating cc per-batch outputs"
validate_files_or_fail "cc batch" "${CC_BATCH_FILES[@]}"

# --- Step 4: hadd into combined files ---
log "hadding bb batches -> ${BB_COMBINED}"
rm -f "${BB_COMBINED}"
# shellcheck disable=SC2086
hadd -f "${BB_COMBINED}" "${BB_DIR}"/single_muon_trees_powheg_bb_fullsim_pp17_part*.root

log "hadding cc batches -> ${CC_COMBINED}"
rm -f "${CC_COMBINED}"
# shellcheck disable=SC2086
hadd -f "${CC_COMBINED}" "${CC_DIR}"/single_muon_trees_powheg_cc_fullsim_pp17_part*.root

log "Validating combined files are valid ROOT files"
validate_files_or_fail "combined" "${BB_COMBINED}" "${CC_COMBINED}"

log "Validating muon_tree is non-empty in bb combined file"
validate_single_muon_tree_nonempty_or_fail "${BB_COMBINED}"

log "Validating muon_tree is non-empty in cc combined file"
validate_single_muon_tree_nonempty_or_fail "${CC_COMBINED}"

# --- Step 5: RDF histogram filling ---
log "Running RDF histogram filling for single-muon (from ${RDF_DIR})"
pushd "${RDF_DIR}" >/dev/null
root -l -b <<'ROOTEOF' || fail "ROOT exited non-zero during single-muon RDF hist filling"
.L RDFBasedHistFillingPowhegFullsimSingleMuon.cxx+
{
    RDFBasedHistFillingPowhegFullsimSingleMuon fs;
    fs.Run();
}
.q
ROOTEOF
popd >/dev/null

# --- Step 6: validate histogram output ---
log "Validating single-muon histogram output: ${SINGLE_MUON_HIST}"
validate_files_or_fail "single-muon histogram" "${SINGLE_MUON_HIST}"

# --- Step 7: plot for medium and tight WPs ---
log "Running plotter for medium + tight WPs (from ${ANALYSIS_DIR})"
pushd "${ANALYSIS_DIR}" >/dev/null
root -l -b <<'ROOTEOF' || fail "ROOT exited non-zero during single-muon plotting"
.L PowhegFullsimDetRespPlotterSingleMuon.cxx+
{
    gROOT->SetBatch(kTRUE);
    PowhegFullsimDetRespPlotterSingleMuon pl_medium(17, false);
    pl_medium.Run();
    PowhegFullsimDetRespPlotterSingleMuon pl_tight(17, true);
    pl_tight.Run();
}
.q
ROOTEOF
popd >/dev/null

log "Single-muon pipeline completed successfully."
