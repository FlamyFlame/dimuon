#!/usr/bin/env bash
set -Eeuo pipefail

# Master PbPb pipeline: runs shared event selection once, then launches
# crossx and trig_eff pipelines in parallel (both skip event selection).
#
# Usage:
#   ./run_pbpb_all.sh
#
# Optional env vars (passed through to sub-pipelines):
#   POLL_SECONDS=45
#   CONDOR_TIMEOUT_SECONDS=0
#   YEARS="23 24 25"
#   SKIP_CONDOR=1              # skip event sel + condor in both pipelines
#   SKIP_EVSEL=1               # skip event selection only
#   RDF_NTHREADS=2             # for trig_eff pipeline

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
EVSEL_DIR="${ANALYSIS_DIR}/plotting_codes/event_selection"

SKIP_CONDOR="${SKIP_CONDOR:-0}"
SKIP_EVSEL="${SKIP_EVSEL:-${SKIP_CONDOR}}"
YEARS=(${YEARS:-23 24 25})
DATA_BASE="/usatlas/u/yuhanguo/usatlasdata/dimuon_data"

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

validate_root_file_quick() {
  [[ -f "$1" ]] && [[ $(stat -c%s "$1" 2>/dev/null || echo 0) -gt 100 ]]
}

# ------ Stage 0: Shared event selection ------
if [[ "$SKIP_EVSEL" -eq 1 ]]; then
  log "SKIP_EVSEL=1 — skipping event selection (reusing existing cuts)"
else
  log "Running shared event selection for all years"
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
  log "Shared event selection complete for all years"
fi

# ------ Launch trig_eff then crossx, SERIALIZED (skip event selection) ------
log "Launching trig_eff then crossx (SERIALIZED -- crossx consumes trig_eff's turn-on fits)"

export SKIP_EVSEL=1
CROSSX_LOG="${SCRIPT_DIR}/crossx_$$.log"
TRIGEFF_LOG="${SCRIPT_DIR}/trigeff_$$.log"

# SERIALIZED, trig-eff FIRST (2026-09-09). These two used to run CONCURRENTLY, which is a race:
# pipeline_pbpb_trig_eff.sh Stage 6 WRITES the single-muon turn-on TF1 fits
#   <year>/trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit.root
# and pipeline_pbpb_crossx.sh Stage 5 READS them to build the per-pair trigger weight. Run in
# parallel, crossx picks up whatever fits happened to be on disk -- i.e. the PREVIOUS round's.
# That is not hypothetical here: the Pb+Pb fit files key their top q*eta bin `_2_00_TO_2_30`
# while the reader now builds `_2_00_TO_2_20` (the coarse edge moved with the fiducial cut), so
# a concurrent run is GUARANTEED to read stale fits and throw inside the RDF event loop -- where
# ROOT swallows the exception and still exits 0, leaving a fresh near-empty output. That failure
# mode already destroyed pbpb_2024/histograms_real_pairs_..._nominal.root (851 bytes, 0 keys).
bash "${SCRIPT_DIR}/pipeline_pbpb_trig_eff.sh" > "$TRIGEFF_LOG" 2>&1
TRIGEFF_RC=$?
log "trig_eff pipeline finished (rc=$TRIGEFF_RC, log: $TRIGEFF_LOG)"

FAIL=0
if [[ $TRIGEFF_RC -ne 0 ]]; then
  log "ERROR: trig_eff pipeline failed (see $TRIGEFF_LOG)"
  log "ABORTING before crossx: it would consume the turn-on fits this stage was supposed to write."
  exit 1
fi
log "trig_eff pipeline completed successfully"

bash "${SCRIPT_DIR}/pipeline_pbpb_crossx.sh" > "$CROSSX_LOG" 2>&1
CROSSX_RC=$?
log "crossx pipeline finished (rc=$CROSSX_RC, log: $CROSSX_LOG)"

if [[ $CROSSX_RC -ne 0 ]]; then
  log "ERROR: crossx pipeline failed (see $CROSSX_LOG)"
  FAIL=1
else
  log "crossx pipeline completed successfully"
fi

if [[ $FAIL -eq 1 ]]; then
  fail "The crossx pipeline failed — check logs above"
fi

log "All PbPb pipelines completed successfully"
