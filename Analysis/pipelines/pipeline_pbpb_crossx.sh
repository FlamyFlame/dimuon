#!/usr/bin/env bash
set -Eeuo pipefail

# End-to-end PbPb crossx/nominal analysis pipeline:
#   0) event selection: derive cuts + produce all event selection plots
#   1) submit NTuple processing condor jobs (nominal mode, all 4 PbPb years)
#   2) wait for all clusters to finish
#   3) validate per-batch ROOT outputs
#   4) hadd per-batch outputs into combined muon_pairs + hists_cut_acceptance per year
#   5) run RDF crossx hist filling per year + validate RDF outputs
#   6) run crossx plotting (combined, auto-discovers all years)
#   7) run before/after trigger efficiency correction sanity check plots
#
# Usage:
#   ./pipeline_pbpb_crossx.sh
#
# Optional env vars:
#   POLL_SECONDS=45
#   CONDOR_TIMEOUT_SECONDS=0   # 0 => no timeout
#   YEARS="23 24 25 26"        # override which years to process
#   SKIP_CONDOR=1              # skip condor submit+wait, reuse existing NTuple outputs
#   SKIP_EVSEL=1               # skip event selection only, still run condor+RDF+plots

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
NTP_DIR="${ANALYSIS_DIR}/NTupleProcessingCode"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
PLOT_DIR="${ANALYSIS_DIR}/plotting_codes/single_b_analysis"
EVSEL_DIR="${ANALYSIS_DIR}/plotting_codes/event_selection"

POLL_SECONDS="${POLL_SECONDS:-45}"
CONDOR_TIMEOUT_SECONDS="${CONDOR_TIMEOUT_SECONDS:-0}"
SKIP_CONDOR="${SKIP_CONDOR:-0}"
SKIP_EVSEL="${SKIP_EVSEL:-${SKIP_CONDOR}}"
YEARS=(${YEARS:-23 24 25 26})
DATA_BASE="/usatlas/u/yuhanguo/usatlasdata/dimuon_data"

# ─── `root -l -b`, NEVER `root -l -b -q`, when the macro comes from a HEREDOC ───────────────
# `root -l -b -q` with no macro argument QUITS BEFORE READING STDIN: the heredoc is never
# executed and the command exits 0 UNCONDITIONALLY, so every check built on it silently PASSES --
# missing trees, zero keys, a zombie file, all of it. The `$(...)` capture form is just as dead:
# it returns an EMPTY string and a 0 status. This trap is documented in run_dr_correction_fits.sh
# (found and fixed there 2026-08-11); the pipelines below still carried it, so their validation
# layer had never actually run. Re-verified 2026-09-09: `root -l -b -q` fed `gSystem->Exit(7)` on
# stdin exits 0 and prints nothing; `root -l -b` exits 7.
# ───────────────────────────────────────────────────────────────────────────────────────────────
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
  root -l -b <<EOF >/dev/null 2>&1
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
  if check_output="$(root -l -b <<EOF
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
# [26]=5 is a PLACEHOLDER: the 2026 skim is submitted as 5 grid tasks
# (SkimCode/run_26hi/InDstxt_PbPb2026_5p36TeV_part1..5.txt), and grid_monitor's chunked
# hadd can split a task into extra part files, so the real count can be larger.  Set it
# from ~/usatlasdata/dimuon_data/pbpb_2026/data_pbpb26_part*.root and keep it EQUAL to
# file_batch_max{26} in NTupleProcessingCode/PbPbExtras.c and to `queue N` in the
# run_pbpb_26*.sub files.  Too few silently processes only part of the 2026 data.
declare -A QUEUE_COUNTS=( [23]=4 [24]=2 [25]=6 [26]=5 )

get_year_dir() { echo "${DATA_BASE}/pbpb_20$1"; }

# NTuple output filename pattern for nominal mode:
#   muon_pairs_pbpb_20YY_partN_single_mu4_mindR_0_02.root
#   hists_cut_acceptance_pbpb_20YY_partN_single_mu4_mindR_0_02.root
get_batch_muon_pairs() {
  local yr="$1" part="$2"
  echo "$(get_year_dir "$yr")/muon_pairs_pbpb_20${yr}_part${part}_single_mu4_mindR_0_02.root"
}
get_batch_hists() {
  local yr="$1" part="$2"
  echo "$(get_year_dir "$yr")/hists_cut_acceptance_pbpb_20${yr}_part${part}_single_mu4_mindR_0_02.root"
}

# Combined (hadd) output filenames
get_combined_muon_pairs() {
  echo "$(get_year_dir "$1")/muon_pairs_pbpb_20${1}_single_mu4_mindR_0_02.root"
}
get_combined_hists() {
  echo "$(get_year_dir "$1")/hists_cut_acceptance_pbpb_20${1}_single_mu4_mindR_0_02.root"
}

# RDF output filename
get_rdf_output() {
  echo "$(get_year_dir "$1")/histograms_real_pairs_pbpb_20${1}_single_mu4_no_trg_plots_nominal.root"
}

# ==============================
# Pipeline execution
# ==============================
log "PbPb crossx/nominal pipeline — years: ${YEARS[*]}"
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
    log "Submitting NTuple nominal jobs for PbPb 20${yr}"
    pushd "$NTP_DIR" >/dev/null
    out="$(condor_submit "run_pbpb_${yr}_nominal.sub")"
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

# ------ Stage 5: RDF crossx hist filling per year + validate ------
# FRESHNESS + a required FILLED histogram, mirroring the pp twin -- not `validate_files_or_fail`
# alone, which only asks for a non-empty file with >=1 key. This is the exact stage that produced
# pbpb_2024/histograms_real_pairs_..._nominal.root at 851 bytes with 0 keys on 2026-09-06: ROOT's
# TRint CATCHES a C++ exception thrown inside the RDF event loop, prints it, and still exits 0, so
# the `|| fail` above cannot see it. A key count also passes on a partially flushed file and on a
# file the RDF never rewrote at all, which is the more insidious case -- it would carry the
# PREVIOUS selection.
for yr in "${YEARS[@]}"; do
  log "Running RDF crossx hist filling for PbPb 20${yr}"
  rdf_stamp="$(mktemp)"
  pushd "$RDF_DIR" >/dev/null
  ./run_crossx_hist_filling_pbpb${yr}.sh || {
    popd >/dev/null
    rm -f "$rdf_stamp"
    fail "RDF crossx hist filling failed for year ${yr}"
  }
  popd >/dev/null
  rdf_out="$(get_rdf_output "$yr")"
  validate_files_or_fail "RDF crossx yr${yr}" "$rdf_out"
  if [[ ! "$rdf_out" -nt "$rdf_stamp" ]]; then
    rm -f "$rdf_stamp"
    fail "RDF crossx output ${rdf_out} is OLDER than this run — the hist filling threw and ROOT swallowed it, so this is the PREVIOUS production's file. Check the log for 'runtime_error'."
  fi
  rm -f "$rdf_stamp"
  # The centrality-inclusive-most bin is present in every Pb+Pb year and is what the combined
  # plotter reads first, so it is the right liveness probe.
  # NO `-q` HERE, and that is not a style choice: `root -l -b -q` with no macro argument QUITS
# BEFORE READING STDIN, so the heredoc is never executed and the command exits 0 UNCONDITIONALLY --
# the probe silently passes, missing histogram and all. That trap is documented and was fixed in
# run_dr_correction_fits.sh on 2026-08-11; these two crossx pipelines still carried it, so this
# check has never actually run. Re-verified 2026-09-09: `root -l -b -q` fed `gSystem->Exit(7)` on
# stdin exits 0 and never prints; `root -l -b` exits 7.
root -l -b >/dev/null 2>&1 <<EOF || fail "RDF crossx output ${rdf_out} has no FILLED h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_ctr0_5 — the event loop threw (ROOT exits 0 anyway)."
TFile* f = TFile::Open("${rdf_out}", "READ");
if (!f || f->IsZombie()) gSystem->Exit(2);
TH1* h = (TH1*)f->Get("h2d_op_crossx_w_signal_cuts_vs_pair_eta_vs_pair_pt_ctr0_5");
if (!h) { f->Close(); gSystem->Exit(3); }
if (h->GetEntries() <= 0) { f->Close(); gSystem->Exit(4); }
f->Close(); gSystem->Exit(0);
EOF
  log "  yr${yr} crossx output is fresh and its signal-region histogram is filled  OK"
done

# ------ Stage 6: Crossx plotting (combined) ------
# The NOMINAL pair-pT view (user decision 2026-09-08) is the UNSUFFIXED histogram family, booked
# on ParamsSet::pT_bins_150 = 16 log bins 9 -> 150 GeV, written to
# pbpb_<years>_combined/{TAA_weighted,counts}/. That is the only fine axis that nests 2:1 inside
# the coarse correction cells (ParamsSet::pair_pt_coarse_bins, 8 log bins over the SAME 9 -> 150
# range) — and, decisively for Pb+Pb, the R_AA global 3Ds belong to that family, so the R_AA input
# and the cross-section share ONE binning. also_pt_120=true additionally refreshes the OPT-IN
# *_combined_pt_120/ variant (the "pt_120" family on ParamsSet::pT_bins_120, 9 -> 120 GeV): an
# alternative VIEW of the same measurement, never a second binning of the nominal result
# (.claude/CLAUDE.md §Binnings). It must be regenerated in the SAME run as the nominal, or the two
# drift apart — they did, 2026-06-19 to 2026-08-04.
log "Running crossx plotting (all years combined; nominal 9-150 GeV + opt-in pt_120 variant)"
pushd "$PLOT_DIR" >/dev/null
root -l -b -q 'plot_single_b_crossx_pbpb.cxx(true)'
popd >/dev/null

# ------ Stage 7: Trigger efficiency correction sanity check ------
log "Running before/after trigger efficiency correction sanity plots"
pushd "$PLOT_DIR" >/dev/null
root -l -b -q 'plot_crossx_trig_corr_sanity.C()'
popd >/dev/null

log "PbPb crossx/nominal pipeline completed successfully for years: ${YEARS[*]}"
