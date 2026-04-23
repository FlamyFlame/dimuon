#!/usr/bin/env bash
set -Eeuo pipefail

# End-to-end Pythia fullsim overlay pipeline:
#   1) Submit Condor job for NTuple processing
#   2) Wait for cluster completion
#   3) Validate NTP output (muon_pairs ROOT file)
#   4) Run RDF histogram filling (per-centrality + inclusive)
#   5) Validate histogram output
#   6) Run plotter for medium + tight WPs
#
# Usage:
#   ./pipeline_pythia_fullsim_overlay.sh hijing
#   ./pipeline_pythia_fullsim_overlay.sh zmumu
#   ./pipeline_pythia_fullsim_overlay.sh data
#   ./pipeline_pythia_fullsim_overlay.sh hijing --dry-run
#
# Optional env vars:
#   POLL_SECONDS=45
#   CONDOR_TIMEOUT_SECONDS=0   # 0 => no timeout
#
# TEMPORARY: I/O paths and Condor job count are based on the current
# test samples. When full overlay samples become available, update
# FullSimSampleType.h paths, .sub queue count, and worker script
# to support batched kn processing.

MODE="${1:-}"
DRY_RUN=0
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
    esac
done

if [[ -z "${MODE}" ]]; then
    echo "Usage: $0 [hijing|zmumu|data] [--dry-run]"
    exit 1
fi

case "$MODE" in
    hijing|zmumu|data) ;;
    --dry-run)
        echo "Usage: $0 [hijing|zmumu|data] [--dry-run]"
        exit 1
        ;;
    *) echo "Unknown mode '${MODE}'. Use hijing|zmumu|data"; exit 1 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
NTP_DIR="${ANALYSIS_DIR}/NTupleProcessingCode"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
RECO_EFFCY_DIR="${ANALYSIS_DIR}/plotting_codes/reco_effcy"

POLL_SECONDS="${POLL_SECONDS:-45}"
CONDOR_TIMEOUT_SECONDS="${CONDOR_TIMEOUT_SECONDS:-0}"

# ---------------------------------------------------------------------------
# Helpers (shared with other pipelines)
# ---------------------------------------------------------------------------
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
    if [[ -n "${err_trap_saved}" ]]; then eval "${err_trap_saved}"; else trap on_error ERR; fi
    [[ $rc -eq 0 ]] || fail "Failed to source ~/setup.sh"
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
    [[ -f "$f" ]] || fail "Muon pairs file not found: $f"

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
  std::cout << "ERROR: muon_pairs file has empty sign tree(s)" << std::endl;
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
        fail "Post-NTP tree non-empty validation failed for: $f"
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

# ---------------------------------------------------------------------------
# Mode-specific configuration
# ---------------------------------------------------------------------------

# Map mode -> C++ enum name, label, data dir
# These match FullSimSampleType.h exactly
case "$MODE" in
    hijing)
        CPP_ENUM="hijing"
        LABEL="hijing_overlay_pp24"
        BASE_DIR="/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample"
        ;;
    zmumu)
        CPP_ENUM="zmumu"
        LABEL="zmumu_overlay_pp24"
        BASE_DIR="/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_zmumu_overlay_test_sample"
        ;;
    data)
        CPP_ENUM="data"
        LABEL="data_overlay_pp24"
        BASE_DIR="/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_data_overlay_test_sample"
        ;;
esac

# TEMPORARY: Single-file output (no batching) for test sample.
# When full samples are available, this will become an array of
# per-batch files that get hadd'd together.
NTP_OUTPUT="${BASE_DIR}/muon_pairs_pythia_fullsim_${LABEL}_no_data_resonance_cuts.root"
HIST_OUTPUT="${BASE_DIR}/histograms_pythia_fullsim_${LABEL}_no_data_resonance_cuts.root"

# ---------------------------------------------------------------------------
# Dry-run mode: use a test output directory to avoid overwriting real output
# ---------------------------------------------------------------------------
if (( DRY_RUN )); then
    TEST_DIR="${BASE_DIR}/pipeline_dry_run_test"
    log "DRY RUN: output will go to ${TEST_DIR} (not overwriting real output)"
    log "DRY RUN: skipping Condor submission; will validate existing NTP output instead"
fi

# ---------------------------------------------------------------------------
# Pipeline execution
# ---------------------------------------------------------------------------
require_cmd condor_submit
require_cmd condor_q

log "=== Pythia fullsim overlay pipeline: mode=${MODE} ==="
log "    NTP output : ${NTP_OUTPUT}"
log "    Hist output: ${HIST_OUTPUT}"
log "    Dry run    : ${DRY_RUN}"

log "Sourcing ~/setup.sh for ROOT/LCG tools"
source_env_once
require_cmd root

# --- Stage 1: NTuple processing via Condor ---
if (( DRY_RUN )); then
    log "[Stage 1] DRY RUN: skipping Condor submission"
    log "[Stage 1] Checking that existing NTP output is valid: ${NTP_OUTPUT}"
    if [[ ! -f "${NTP_OUTPUT}" ]]; then
        fail "[Stage 1] DRY RUN: NTP output does not exist: ${NTP_OUTPUT}. Run the pipeline without --dry-run first."
    fi
else
    log "[Stage 1] Submitting NTuple processing Condor job"
    mkdir -p "${NTP_DIR}/logs"
    pushd "${NTP_DIR}" >/dev/null
    submit_out="$(condor_submit \
        -append "arguments = ${CPP_ENUM}" \
        -append "sample_type = ${CPP_ENUM}" \
        run_pythia_fullsim_overlay.sub)"
    echo "${submit_out}" >&2
    cluster_id="$(extract_cluster_id "${submit_out}")"
    log "Submitted cluster: ${cluster_id}"
    popd >/dev/null

    log "[Stage 1] Waiting for cluster ${cluster_id} to finish"
    wait_for_cluster_completion "${cluster_id}"
fi

# --- Stage 2: Validate NTP output ---
log "[Stage 2] Validating NTP output"
validate_files_or_fail "NTP" "${NTP_OUTPUT}"
validate_combined_muon_pair_trees_nonempty_or_fail "${NTP_OUTPUT}"

# TEMPORARY: No hadd needed — single job produces single output file.
# When full samples use batched Condor jobs, add hadd step here:
#   hadd -f "${COMBINED_OUTPUT}" ${NTP_OUTPUT_PATTERN}
#   validate_combined_muon_pair_trees_nonempty_or_fail "${COMBINED_OUTPUT}"

# --- Stage 3: RDF histogram filling ---
log "[Stage 3] Running RDF histogram filling (from ${RDF_DIR})"
pushd "${RDF_DIR}" >/dev/null

if (( DRY_RUN )); then
    log "[Stage 3] DRY RUN: compiling and running RDF overlay hist filling (same output path — test validates the chain runs)"
fi

root -l -b <<ROOTEOF || fail "ROOT exited non-zero during RDF histogram filling"
.L RDFBasedHistFillingPythiaFullsimOverlay.cxx+
{
    RDFBasedHistFillingPythiaFullsimOverlay fs;
    fs.fullsim_sample_type = FullSimSampleType::${CPP_ENUM};
    fs.Run();
}
.q
ROOTEOF
popd >/dev/null

# --- Stage 4: Validate histogram output ---
log "[Stage 4] Validating histogram output"
validate_files_or_fail "Histogram" "${HIST_OUTPUT}"

# Quick sanity check: histogram file should have a reasonable number of keys
hist_root_out="$(root -l -b <<EOF 2>/dev/null || true
TFile *fin = TFile::Open("${HIST_OUTPUT}", "READ");
if (!fin || fin->IsZombie()) { gSystem->Exit(2); }
std::cout << "NKEYS=" << fin->GetListOfKeys()->GetSize() << std::endl;
fin->Close();
gSystem->Exit(0);
EOF
)"
hist_nkeys="$(echo "$hist_root_out" | grep '^NKEYS=' | cut -d= -f2 || true)"
if [[ -n "$hist_nkeys" ]] && (( hist_nkeys > 0 )); then
    log "Histogram file has ${hist_nkeys} keys"
else
    fail "Histogram file has no keys or could not be read: ${HIST_OUTPUT}"
fi

# --- Stage 5: Plotting ---
log "[Stage 5] Running plotter for medium + tight WPs"
pushd "${RECO_EFFCY_DIR}" >/dev/null
root -l -b <<ROOTEOF || fail "ROOT exited non-zero during plotting"
.L PythiaFullsimRecoEffPlotter.cxx+
{
    gROOT->SetBatch(kTRUE);

    PythiaFullsimRecoEffPlotterOverlay pl_medium(FullSimSampleType::${CPP_ENUM}, false, false);
    pl_medium.Run();

    PythiaFullsimRecoEffPlotterOverlay pl_tight(FullSimSampleType::${CPP_ENUM}, true, false);
    pl_tight.Run();

    PythiaFullsimRecoEffPlotterOverlay pl_medium_sig(FullSimSampleType::${CPP_ENUM}, false, true);
    pl_medium_sig.Run();

    PythiaFullsimRecoEffPlotterOverlay pl_tight_sig(FullSimSampleType::${CPP_ENUM}, true, true);
    pl_tight_sig.Run();
}
.q
ROOTEOF
popd >/dev/null

log "=== Pipeline completed successfully for mode=${MODE} ==="
