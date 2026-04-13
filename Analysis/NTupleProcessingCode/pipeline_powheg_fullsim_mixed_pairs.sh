#!/usr/bin/env bash
# End-to-end mixed-pair efficiency pipeline:
#
#   0) [Optional] smoke_test - run a 200-pair bb mixing test before the full job
#   1) Submit Condor mixing jobs (480 batches, bb+cc combined input)
#   2) Wait for the cluster to finish
#   3) Validate all 480 batch output files (existence + non-zero size);
#      spot-check tree non-emptiness on a representative subset
#   4) Run RDF histogram filling for mixed muon pairs
#   5) Validate the histogram output file
#   6) Run the plotter for both medium and tight working points
#
# Prerequisite: pipeline_powheg_fullsim_single_muon.sh has been run and
#   the corrected combined single-muon trees exist at:
#     .../bb.../single_muon_trees_powheg_bb_fullsim_pp17.root
#     .../cc.../single_muon_trees_powheg_cc_fullsim_pp17.root
#
# Usage:
#   ./pipeline_powheg_fullsim_mixed_pairs.sh [--smoke-test] [--mass-filter|--no-mass-filter] [--skip-mixing]
#
#   --smoke-test      : run the 200-pair bb smoke test before the full pipeline
#   --mass-filter     : apply truth_minv in [1,3] GeV filter (default); files go to
#                       mixed_mass_1_3GeV/, histograms named *_mixed_mass_1_3GeV.root,
#                       plots in run2_reco_effcy_plots_mass_1_3GeV/
#   --no-mass-filter  : no minv filter; files go to mixed/, histograms named
#                       *_mixed.root, plots in run2_reco_effcy_plots/
#   --skip-mixing     : skip Condor job submission (Steps 1-2); go straight to
#                       validating existing mixed files then RDF filling + plotting.
#                       Useful when mixing jobs are already done.
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
# How long (seconds) a lone remaining job may run after all others finished
# before it is killed and skipped.  0 = no timeout.  Default: 12 h.
LONE_JOB_TIMEOUT_SECONDS="${LONE_JOB_TIMEOUT_SECONDS:-43200}"
RUN_SMOKE_TEST=0
# 1 = apply truth_minv in [1,3] GeV filter (default); 0 = no filter
APPLY_MASS_FILTER=1
# 1 = skip Condor job submission (Steps 1-2) and go straight to validation (Step 3)
SKIP_MIXING=0

# Batch indices (1-based) killed/skipped during waiting; populated by
# wait_for_cluster_completion and consumed by later validation and RDF steps.
SKIPPED_BATCHES=()

for arg in "$@"; do
    case "$arg" in
        --smoke-test)      RUN_SMOKE_TEST=1 ;;
        --mass-filter)     APPLY_MASS_FILTER=1 ;;
        --no-mass-filter)  APPLY_MASS_FILTER=0 ;;
        --skip-mixing)     SKIP_MIXING=1 ;;
        *) echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Paths — derived from APPLY_MASS_FILTER
# ---------------------------------------------------------------------------
POWHEG_BASE="/usatlas/u/yuhanguo/usatlasdata/powheg_full_sample"
BB_DIR="${POWHEG_BASE}/user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp17.bb.Feb2026.v1._MYSTREAM"
CC_DIR="${POWHEG_BASE}/user.yuhang.TrigRates.dimuon.PowhegPythia.fullsim.pp17.cc.Feb2026.v1._MYSTREAM"

if (( APPLY_MASS_FILTER )); then
    MIXED_SUBDIR="mixed_mass_1_3GeV"
    MIXED_HIST_SUFFIX="_mass_1_3GeV"
    CPP_APPLY_MASS_FILTER="true"
else
    MIXED_SUBDIR="mixed"
    MIXED_HIST_SUFFIX=""
    CPP_APPLY_MASS_FILTER="false"
fi

MIXED_DIR="${POWHEG_BASE}/${MIXED_SUBDIR}"

BB_COMBINED="${BB_DIR}/single_muon_trees_powheg_bb_fullsim_pp17.root"
CC_COMBINED="${CC_DIR}/single_muon_trees_powheg_cc_fullsim_pp17.root"

MIXED_HIST="${POWHEG_BASE}/histograms_powheg_fullsim_pp17_mixed${MIXED_HIST_SUFFIX}.root"

# Spot-check these batch numbers for full tree validation (first / mid / last)
SPOT_CHECK_BATCHES=(1 120 240 360 480)

# ---------------------------------------------------------------------------
# Helpers
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
    local t0 lone_since lone_proc_id
    t0="$(date +%s)"
    lone_since=0
    lone_proc_id=-1

    while true; do
        local qout
        qout="$(condor_q "${cluster_id}" -autoformat ClusterId ProcId JobStatus 2>/dev/null || true)"

        if [[ -z "$qout" ]]; then
            log "Cluster ${cluster_id} no longer in condor_q (finished/left queue)."
            break
        fi

        local total=0 idle=0 running=0 held=0 other=0 cur_running_proc=-1
        while read -r cid pid st; do
            [[ -n "${cid:-}" ]] || continue
            total=$((total + 1))
            case "$st" in
                1) idle=$((idle + 1)) ;;
                2) running=$((running + 1)); cur_running_proc="${pid}" ;;
                5) held=$((held + 1)) ;;
                *) other=$((other + 1)) ;;
            esac
        done <<< "$qout"

        log "Cluster ${cluster_id}: queued=${total}, idle=${idle}, running=${running}, held=${held}, other=${other}"

        if (( held > 0 )); then
            fail "Cluster ${cluster_id} has held jobs; inspect with: condor_q ${cluster_id} -hold"
        fi

        # ---- Lone-job tracking ----
        # When exactly one job remains (and none are still idle/pending),
        # start a per-lone-job timer.  If it exceeds LONE_JOB_TIMEOUT_SECONDS,
        # check whether its output file already exists (job completed its work
        # but ROOT hung during teardown); if so, kill it and skip it
        # downstream; otherwise kill it and fail.
        if (( running == 1 && idle == 0 && LONE_JOB_TIMEOUT_SECONDS > 0 )); then
            if (( lone_since == 0 )); then
                lone_since="$(date +%s)"
                lone_proc_id="${cur_running_proc}"
                log "Lone-job detected: ${cluster_id}.${lone_proc_id} -- lone-job timer started (timeout=${LONE_JOB_TIMEOUT_SECONDS}s)"
            fi
            local lone_elapsed=$(( $(date +%s) - lone_since ))
            if (( lone_elapsed > LONE_JOB_TIMEOUT_SECONDS )); then
                local lone_batch=$(( lone_proc_id + 1 ))
                local lone_output="${MIXED_DIR}/muon_pairs_powheg_bbcc_fullsim_mixed_batch${lone_batch}.root"
                log "Lone-job timeout: ${cluster_id}.${lone_proc_id} (batch ${lone_batch}) has been the sole remaining job for ${lone_elapsed}s (>${LONE_JOB_TIMEOUT_SECONDS}s)"
                if [[ -f "${lone_output}" && -s "${lone_output}" ]]; then
                    local fsize
                    fsize="$(du -h "${lone_output}" | cut -f1)"
                    log "Output file exists (${fsize}): ${lone_output}. Job finished its work but ROOT hung during teardown. Killing and treating batch ${lone_batch} as complete (skipped in validation + RDF)."
                    condor_rm "${cluster_id}.${lone_proc_id}" || true
                    SKIPPED_BATCHES+=("${lone_batch}")
                    break
                else
                    log "Output file missing or empty: ${lone_output}. Killing lone job ${cluster_id}.${lone_proc_id} (batch ${lone_batch}) -- it produced no output."
                    condor_rm "${cluster_id}.${lone_proc_id}" || true
                    fail "Lone job ${cluster_id}.${lone_proc_id} (batch ${lone_batch}) timed out with no output."
                fi
            fi
        else
            # Reset lone tracker if more jobs become active (unlikely but safe)
            if (( running != 1 || idle != 0 )); then
                lone_since=0
                lone_proc_id=-1
            fi
        fi

        # ---- Global cluster timeout ----
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

    if (( ${#SKIPPED_BATCHES[@]} > 0 )); then
        log "Batches killed/skipped due to lone-job timeout: ${SKIPPED_BATCHES[*]}"
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

validate_mixed_pair_trees_nonempty_or_fail() {
    local f="$1"
    [[ -f "$f" ]] || fail "Mixed pair file not found: $f"

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
  std::cout << "ERROR: missing sign tree(s) in $f" << std::endl;
  fin->Close();
  gSystem->Exit(3);
}
Long64_t n1 = t1->GetEntries();
Long64_t n2 = t2->GetEntries();
std::cout << "sign1=" << n1 << " sign2=" << n2 << " in $f" << std::endl;
if (n1 <= 0 || n2 <= 0) {
  std::cout << "ERROR: empty sign tree(s) in $f" << std::endl;
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
        fail "Mixed pair tree non-empty validation failed for: $f"
    fi
}

# ---------------------------------------------------------------------------
# Smoke test: mix 200 bb pairs from the corrected single-muon trees and
# confirm the output file has non-empty sign trees.  Exits on failure so
# we never submit 480 condor jobs with broken input.
# ---------------------------------------------------------------------------
run_smoke_test() {
    log "[SMOKE TEST] Starting 200-pair bb mixing smoke test (apply_mass_filter=${APPLY_MASS_FILTER})"

    local smoke_output="${MIXED_DIR}/muon_pairs_powheg_bbcc_fullsim_mixed_batch999_smoke.root"
    mkdir -p "${MIXED_DIR}"

    pushd "${SCRIPT_DIR}" >/dev/null
    root -l -b <<ROOTEOF || fail "[SMOKE TEST] ROOT exited non-zero during smoke mixing"
.L mix_powheg_single_muon_pairs.C
run_powheg_single_muon_pair_mixing(
    "bb",       // mc_mode
    200,        // target pairs
    999,        // file_batch
    17,         // run_year
    "",         // input_file = "" => use AddDefaultInputFiles (bb+cc combined)
    "${smoke_output}",
    999,        // rng_seed
    ${CPP_APPLY_MASS_FILTER}  // apply_mass_filter
);
.q
ROOTEOF
    popd >/dev/null

    log "[SMOKE TEST] Validating smoke output: ${smoke_output}"
    validate_mixed_pair_trees_nonempty_or_fail "${smoke_output}"

    log "[SMOKE TEST] PASSED -- removing smoke output file"
    rm -f "${smoke_output}"
}

# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------
require_cmd condor_submit
require_cmd condor_q

log "Sourcing ~/setup.sh for ROOT / condor tools"
source_env_once
require_cmd root

log "=== Pipeline mode: APPLY_MASS_FILTER=${APPLY_MASS_FILTER}  SKIP_MIXING=${SKIP_MIXING} ==="
log "    mixed subdir : ${MIXED_SUBDIR}"
log "    histogram    : ${MIXED_HIST}"
log "    plot prefix  : run2_reco_effcy_plots${MIXED_HIST_SUFFIX}/"

if (( SKIP_MIXING )); then
    log "=== --skip-mixing: skipping Condor job submission and waiting (Steps 1-2) ==="
    log "    Will validate existing files in ${MIXED_DIR}"
else
    # --- Prerequisite check: corrected single-muon combined trees must exist ---
    log "Checking prerequisite: combined single-muon trees from pipeline_powheg_fullsim_single_muon.sh"
    [[ -f "${BB_COMBINED}" ]] || fail "bb combined single-muon tree not found: ${BB_COMBINED}. Run pipeline_powheg_fullsim_single_muon.sh first."
    [[ -f "${CC_COMBINED}" ]] || fail "cc combined single-muon tree not found: ${CC_COMBINED}. Run pipeline_powheg_fullsim_single_muon.sh first."
    log "Prerequisites OK."

    # --- Optional smoke test ---
    if (( RUN_SMOKE_TEST )); then
        run_smoke_test
    fi

    # --- Step 1: submit mixing condor cluster ---
    log "Submitting mixing Condor jobs (run_powheg_fullsim_single_muon_mixing.sub, 480 batches, APPLY_MASS_FILTER=${APPLY_MASS_FILTER})"
    mkdir -p "${MIXED_DIR}"
    pushd "${SCRIPT_DIR}" >/dev/null
    # Export APPLY_MASS_FILTER so the Condor wrapper script picks it up.
    mix_submit_out="$(APPLY_MASS_FILTER=${APPLY_MASS_FILTER} condor_submit \
        -append "environment = APPLY_MASS_FILTER=${APPLY_MASS_FILTER}" \
        run_powheg_fullsim_single_muon_mixing.sub)"
    echo "${mix_submit_out}" >&2
    mix_cluster="$(extract_cluster_id "${mix_submit_out}")"
    log "Submitted mixing cluster: ${mix_cluster}"
    popd >/dev/null

    # --- Step 2: wait for cluster ---
    log "Waiting for mixing cluster ${mix_cluster} to finish (480 jobs)"
    wait_for_cluster_completion "${mix_cluster}"
fi

# --- Step 3: validate all 480 batch output files ---
is_skipped_batch() {
    local b="$1"
    local s
    for s in "${SKIPPED_BATCHES[@]}"; do
        [[ "$s" == "$b" ]] && return 0
    done
    return 1
}

if (( ${#SKIPPED_BATCHES[@]} > 0 )); then
    _skipped_str="${SKIPPED_BATCHES[*]}"
else
    _skipped_str="none"
fi
log "Checking existence and non-zero size of all 480 mixed batch files (skipped: ${_skipped_str})"
unset _skipped_str
n_bad=0
for batch in $(seq 1 480); do
    if is_skipped_batch "${batch}"; then
        log "Skipping validation for batch ${batch} (lone-job skip list)"
        continue
    fi
    f="${MIXED_DIR}/muon_pairs_powheg_bbcc_fullsim_mixed_batch${batch}.root"
    if [[ ! -f "$f" ]]; then
        echo "[$(now)] MISSING: $f" >&2
        n_bad=$((n_bad + 1))
    elif [[ ! -s "$f" ]]; then
        echo "[$(now)] EMPTY: $f" >&2
        n_bad=$((n_bad + 1))
    fi
done
[[ $n_bad -eq 0 ]] || fail "${n_bad} mixed batch file(s) missing or empty (excluding ${#SKIPPED_BATCHES[@]} skipped)"
log "All non-skipped mixed batch files exist and are non-empty."

log "Spot-checking tree non-emptiness for batches: ${SPOT_CHECK_BATCHES[*]}"
for batch in "${SPOT_CHECK_BATCHES[@]}"; do
    if is_skipped_batch "${batch}"; then
        log "Skipping spot-check for batch ${batch} (lone-job skip list)"
        continue
    fi
    f="${MIXED_DIR}/muon_pairs_powheg_bbcc_fullsim_mixed_batch${batch}.root"
    validate_mixed_pair_trees_nonempty_or_fail "$f"
    log "Spot-check OK: batch ${batch}"
done

# --- Step 4: RDF histogram filling for mixed pairs ---
# Build C++ snippet to insert skipped batches into fs.skip_batches_mixed
skip_cpp_lines=""
if (( ${#SKIPPED_BATCHES[@]} > 0 )); then
    log "Passing ${#SKIPPED_BATCHES[@]} skipped batch(es) to RDF filler: ${SKIPPED_BATCHES[*]}"
    for _b in "${SKIPPED_BATCHES[@]}"; do
        skip_cpp_lines+="    fs.skip_batches_mixed.insert(${_b});"$'\n'
    done
fi

log "Running RDF histogram filling for mixed muon pairs (from ${RDF_DIR}, mixed_subdir=${MIXED_SUBDIR})"
pushd "${RDF_DIR}" >/dev/null
root -l -b <<ROOTEOF || fail "ROOT exited non-zero during mixed-pair RDF hist filling"
.L RDFBasedHistFillingPowhegFullsim.cxx+
{
    RDFBasedHistFillingPowhegFullsim fs(17, true);  // run_year=17, useMixed=true
    fs.mixed_subdir = "${MIXED_SUBDIR}";
${skip_cpp_lines}    fs.Run();
}
.q
ROOTEOF
popd >/dev/null

# --- Step 5: validate histogram output ---
log "Validating mixed-pair histogram output: ${MIXED_HIST}"
[[ -f "${MIXED_HIST}" && -s "${MIXED_HIST}" ]] || fail "Histogram output not found or empty: ${MIXED_HIST}"
root -l -b -q <<EOF >/dev/null 2>&1
TFile *fin = TFile::Open("${MIXED_HIST}", "READ");
if (!fin || fin->IsZombie()) { gSystem->Exit(2); }
if (!fin->GetListOfKeys() || fin->GetListOfKeys()->GetSize() <= 0) { fin->Close(); gSystem->Exit(3); }
fin->Close();
gSystem->Exit(0);
EOF
log "Histogram output OK: ${MIXED_HIST}"

# --- Step 6: plot for medium and tight WPs ---
log "Running pair plotter for medium + tight WPs (from ${ANALYSIS_DIR}, mixed_hist_suffix=${MIXED_HIST_SUFFIX})"
pushd "${ANALYSIS_DIR}" >/dev/null
root -l -b <<ROOTEOF || fail "ROOT exited non-zero during mixed-pair plotting"
.L RecoEffyRetRespPlotter.cxx+
{
    gROOT->SetBatch(kTRUE);
    RecoEffyRetRespPlotter pl_medium(17, false, false);
    pl_medium.mixed_hist_name_suffix = "${MIXED_HIST_SUFFIX}";
    pl_medium.Run();
    RecoEffyRetRespPlotter pl_tight(17, true, false);
    pl_tight.mixed_hist_name_suffix = "${MIXED_HIST_SUFFIX}";
    pl_tight.Run();
}
.q
ROOTEOF
popd >/dev/null

log "Mixed-pairs pipeline completed successfully."
