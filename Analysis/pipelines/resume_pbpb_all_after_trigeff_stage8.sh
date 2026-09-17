#!/usr/bin/env bash
# Resume run_pbpb_all.sh after a failure in pipeline_pbpb_trig_eff.sh Stage 8 (plotting):
# event selection is DONE and the trig-eff NTuple Condor outputs are on disk, so re-enter the
# trig-eff pipeline with SKIP_EVSEL=1 SKIP_CONDOR=1 (re-hadd + RDF + fits + plots, deterministic),
# then the medium-WP refits, then the crossx pipeline WITH its own Condor stage (SKIP_EVSEL only).
# Same serialization and rc-capture idiom as run_pbpb_all.sh (`cmd || RC=$?`).
set -Eeuo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export YEARS="${YEARS:-23 25 26}"
now() { date '+%F %T'; }
log() { echo "[$(now)] $*"; }

TRIGEFF_LOG="${SCRIPT_DIR}/trigeff_resume_$$.log"
RC=0
SKIP_EVSEL=1 SKIP_CONDOR=1 bash "${SCRIPT_DIR}/pipeline_pbpb_trig_eff.sh" > "$TRIGEFF_LOG" 2>&1 || RC=$?
log "trig_eff pipeline finished (rc=$RC, log: $TRIGEFF_LOG)"
[[ $RC -eq 0 ]] || { log "ERROR: trig_eff failed -- ABORTING before crossx"; exit 1; }

MEDIUM_LOG="${SCRIPT_DIR}/trigeff_medium_resume_$$.log"
RC=0
SAMPLES="${YEARS}" bash "${SCRIPT_DIR}/run_data_trigeff_medium_wp.sh" > "$MEDIUM_LOG" 2>&1 || RC=$?
log "medium-WP trig_eff finished (rc=$RC, log: $MEDIUM_LOG)"
[[ $RC -eq 0 ]] || { log "ERROR: medium-WP trig_eff failed -- ABORTING before crossx"; exit 1; }

CROSSX_LOG="${SCRIPT_DIR}/crossx_resume_$$.log"
RC=0
SKIP_EVSEL=1 SKIP_CONDOR=0 bash "${SCRIPT_DIR}/pipeline_pbpb_crossx.sh" > "$CROSSX_LOG" 2>&1 || RC=$?
log "crossx pipeline finished (rc=$RC, log: $CROSSX_LOG)"
[[ $RC -eq 0 ]] || { log "ERROR: crossx pipeline failed (see $CROSSX_LOG)"; exit 1; }
log "All PbPb pipelines completed successfully"
