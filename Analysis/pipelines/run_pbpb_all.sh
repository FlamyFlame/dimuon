#!/usr/bin/env bash
set -Eeuo pipefail

# Master PbPb pipeline: runs shared event selection once, then runs, IN ORDER,
# trig_eff -> medium-WP turn-on fits -> crossx (both sub-pipelines skip event selection).
# The order is a data dependency, not a preference: crossx consumes the turn-on fits.
#
# Usage:
#   ./run_pbpb_all.sh
#
# Optional env vars (passed through to sub-pipelines):
#   POLL_SECONDS=45
#   CONDOR_TIMEOUT_SECONDS=0
#   YEARS="23 24 25 26"
#   SKIP_CONDOR=1              # skip event sel + condor in both pipelines
#   SKIP_EVSEL=1               # skip event selection only
#   RDF_NTHREADS=2             # for trig_eff pipeline

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
EVSEL_DIR="${ANALYSIS_DIR}/plotting_codes/event_selection"

SKIP_CONDOR="${SKIP_CONDOR:-0}"
SKIP_EVSEL="${SKIP_EVSEL:-${SKIP_CONDOR}}"
# YEARS is consumed TWICE: locally as an array, and by the sub-pipelines as an inherited scalar.
# `YEARS=(${YEARS:-...})` did only the first -- it replaced the scalar with an array, and arrays
# are not exported, so every child saw YEARS UNSET and silently fell back to all three years.
# Consequence during a rerun: `YEARS=24 SKIP_CONDOR=1 ./run_pbpb_all.sh`, the natural way to
# re-enter after one year fails, resubmits ~24 Condor jobs and overwrites 2023 and 2025 outputs
# that were already good. Keep the scalar (exported) and derive the array from it.
YEARS_STR="${YEARS:-23 24 25 26}"
DATA_BASE="/usatlas/u/yuhanguo/usatlasdata/dimuon_data"

# ─── Drop years whose raw skim NTUPs are not on disk yet ───────────────────────────────
# The default year list names every Pb+Pb data-taking period, including one whose grid
# skim may still be running.  Without this filter the pipeline would abort on the first
# stage that touches the absent year -- taking the years that ARE ready down with it.
# A skipped year is announced loudly, once, so it can never be mistaken for a year that
# was processed.  Pass YEARS explicitly to override the discovery.
filter_years_with_data() {
  local kept=() yr
  for yr in "$@"; do
    if compgen -G "${DATA_BASE}/pbpb_20${yr}/data_pbpb${yr}_part*.root" > /dev/null; then
      kept+=("${yr}")
    else
      echo "[SKIP] Pb+Pb 20${yr}: no raw skim NTUPs in ${DATA_BASE}/pbpb_20${yr}/ -- year skipped." >&2
    fi
  done
  # NOTE: this function is called inside $(...) / <(...), i.e. in a SUBSHELL, so it must
  # NOT try to abort the pipeline itself -- an `exit` here would only end the subshell and
  # leave the caller with an EMPTY year list that silently processes nothing.  The caller
  # checks for empty and aborts.
  # The emptiness guard matters: `printf '%s\n'` with ZERO arguments still prints one
  # blank line, which mapfile would turn into a one-element array containing "", and the
  # caller's `${#YEARS[@]} -eq 0` test would then pass with nothing to process.
  if [[ ${#kept[@]} -gt 0 ]]; then
    printf '%s\n' "${kept[@]}"
  fi
}
# Resolve the year list ONCE here and export the result, so every child script
# (and the medium-WP driver below) sees the same, already-filtered set.
mapfile -t YEARS_ARR < <(filter_years_with_data ${YEARS_STR})
if [[ ${#YEARS_ARR[@]} -eq 0 ]]; then
  echo "[FATAL] none of the requested Pb+Pb years has raw skim NTUPs on disk -- nothing to do." >&2
  exit 1
fi
YEARS_STR="${YEARS_ARR[*]}"
export YEARS="${YEARS_STR}"
echo "[INFO] Pb+Pb years to process: ${YEARS_STR}"

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
  for yr in "${YEARS_ARR[@]}"; do
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
# `$?` after a BARE command is unreachable here: line 2 is `set -Eeuo pipefail`, so a non-zero
# child aborts this script before the capture runs, and the operator gets the generic ERR-trap
# line instead of the specific message. The serialization itself was always correct -- crossx is
# a later command and cannot start early -- but its diagnostics were dead code.
# `cmd || RC=$?` with RC preset to 0 is the form that works. NOT `if ! cmd; then RC=$?; fi`:
# inside that `then` branch `$?` is the status of the `!` NEGATION, which is always 0, so the
# real exit code is lost and every failure reports rc=0. Verified both forms empirically before
# committing -- the negation form silently captured 0 from a child that exited 7.
TRIGEFF_RC=0
bash "${SCRIPT_DIR}/pipeline_pbpb_trig_eff.sh" > "$TRIGEFF_LOG" 2>&1 || TRIGEFF_RC=$?
log "trig_eff pipeline finished (rc=$TRIGEFF_RC, log: $TRIGEFF_LOG)"

if [[ $TRIGEFF_RC -ne 0 ]]; then
  log "ERROR: trig_eff pipeline failed (see $TRIGEFF_LOG)"
  log "ABORTING before crossx: it would consume the turn-on fits this stage was supposed to write."
  exit 1
fi
log "trig_eff pipeline completed successfully"

# ------ MEDIUM-WP turn-on fits, before crossx ------
# pipeline_pbpb_trig_eff.sh Stage 5 leaves `isTight` at its default true, so it refreshes only the
# TIGHT tag-and-probe file and its fit. Without this stage the Pb+Pb Medium turn-on fits would be
# the ONLY correction left on the superseded (-1.30,-1.05) gap window and the old `_2_00_TO_2_30`
# q*eta key, while everything around them moved -- and the WP registry requires the Medium leg to
# exist for the WP systematic. It runs BEFORE crossx so a Medium crossx pass can never read a
# Tight-vintage fit. Note this CREATES the 2024/2025 files: neither year had one on disk.
MEDIUM_LOG="${SCRIPT_DIR}/trigeff_medium_$$.log"
MEDIUM_RC=0
SAMPLES="${YEARS_STR}" bash "${SCRIPT_DIR}/run_data_trigeff_medium_wp.sh" > "$MEDIUM_LOG" 2>&1 || MEDIUM_RC=$?
log "medium-WP trig_eff finished (rc=$MEDIUM_RC, log: $MEDIUM_LOG)"
if [[ $MEDIUM_RC -ne 0 ]]; then
  log "ERROR: medium-WP trig_eff failed (see $MEDIUM_LOG)"
  log "ABORTING before crossx: the Medium turn-on fits would stay on the superseded selection."
  exit 1
fi
log "medium-WP trig_eff completed successfully"

CROSSX_RC=0
bash "${SCRIPT_DIR}/pipeline_pbpb_crossx.sh" > "$CROSSX_LOG" 2>&1 || CROSSX_RC=$?
log "crossx pipeline finished (rc=$CROSSX_RC, log: $CROSSX_LOG)"

if [[ $CROSSX_RC -ne 0 ]]; then
  log "ERROR: crossx pipeline failed (see $CROSSX_LOG)"
  fail "The crossx pipeline failed — check $CROSSX_LOG"
fi
log "crossx pipeline completed successfully"

log "All PbPb pipelines completed successfully"
