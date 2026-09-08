#!/usr/bin/env bash
set -Eeuo pipefail

# End-to-end PP crossx/nominal analysis pipeline:
#   1) submit NTuple processing condor jobs (nominal mode, trigger_mode=3 / 2mu4)
#   2) wait for all clusters to finish
#   3) validate per-batch ROOT outputs
#   4) hadd per-batch outputs into combined muon_pairs + hists_cut_acceptance
#   5) run RDF crossx hist filling (trigger_mode=3, nominal mode)
#   6) run crossx plotting
#   7) run before/after trigger efficiency correction sanity check plots
#   8) [optional] run MC-data (POWHEG/Pythia vs pp) comparison plots
#
# PP has NO event selection (no ZDC, no centrality, no nTrk_HITight).
#
# Usage:
#   ./pipeline_pp_crossx.sh
#
# Optional env vars:
#   POLL_SECONDS=45
#   CONDOR_TIMEOUT_SECONDS=0   # 0 => no timeout
#   SKIP_CONDOR=1              # skip condor submit+wait, reuse existing NTuple outputs
#   SKIP_MC_DATA_COMPR=1       # skip MC-data comparison plots (default: 0, i.e. enabled)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
NTP_DIR="${ANALYSIS_DIR}/NTupleProcessingCode"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
PLOT_DIR="${ANALYSIS_DIR}/plotting_codes/single_b_analysis"

POLL_SECONDS="${POLL_SECONDS:-45}"
CONDOR_TIMEOUT_SECONDS="${CONDOR_TIMEOUT_SECONDS:-0}"
SKIP_CONDOR="${SKIP_CONDOR:-0}"
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

# NTuple output filename pattern for PP nominal (trigger_mode=3, 2mu4):
#   muon_pairs_pp_2024_partN_2mu4_mindR_0_02.root
#   hists_cut_acceptance_pp_2024_partN_2mu4_mindR_0_02.root
# The mindR suffix is added at runtime when the skim has mindR branches.
TRIG_SUFFIX="_2mu4_mindR_0_02"

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
  echo "${PP_DIR}/histograms_real_pairs_pp_2024_2mu4_nominal.root"
}

# ==============================
# Pipeline execution
# ==============================
log "PP crossx/nominal pipeline"
log "Sourcing ~/setup.sh"
source_env_once
require_cmd root
require_cmd hadd

if [[ "$SKIP_CONDOR" -eq 1 ]]; then
  log "SKIP_CONDOR=1 — skipping condor submit+wait, reusing existing NTuple outputs"
else
  require_cmd condor_submit
  require_cmd condor_q

  # ------ Stage 1: Submit NTuple nominal condor jobs ------
  mkdir -p "${NTP_DIR}/logs"
  log "Submitting NTuple nominal jobs for PP 2024 (trigger_mode=3, 2mu4)"
  pushd "$NTP_DIR" >/dev/null
  out="$(condor_submit run_pp_24_nominal.sub)"
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

# ------ Stage 5: RDF crossx hist filling + validate ------
# FRESHNESS, not just existence. ROOT's TRint CATCHES a C++ exception thrown inside the event
# loop, prints it, and still exits 0 -- so `run_crossx_hist_filling_pp24.sh` can "succeed" while
# leaving the PREVIOUS run's histogram file untouched on disk. That happened on 2026-08-18 (a
# muon at exactly q*eta = 2.30 had no fitted turn-on): the file opened, had keys, passed
# validation, and the pipeline only failed two stages later in the plotter. Stamp the time before
# the run and require the output to be newer than it.
log "Running RDF crossx hist filling for PP 2024"
rdf_stamp="$(mktemp)"
pushd "$RDF_DIR" >/dev/null
./run_crossx_hist_filling_pp24.sh || {
  popd >/dev/null
  rm -f "$rdf_stamp"
  fail "RDF crossx hist filling failed for PP 2024"
}
popd >/dev/null
rdf_out="$(get_rdf_output)"
validate_files_or_fail "RDF crossx pp24" "$rdf_out"
if [[ ! "$rdf_out" -nt "$rdf_stamp" ]]; then
  rm -f "$rdf_stamp"
  fail "RDF crossx output ${rdf_out} is OLDER than this run — the hist filling threw and ROOT swallowed it. Check the log for 'runtime_error'."
fi
rm -f "$rdf_stamp"
# ...and the file must actually CONTAIN the crossx spectrum. A throw mid-event-loop leaves a
# freshly-RECREATEd near-empty file, which is newer than the stamp and still opens.
root -l -b -q <<EOF >/dev/null 2>&1
TFile* f = TFile::Open("${rdf_out}", "READ");
if (!f || f->IsZombie()) gSystem->Exit(2);
if (!f->Get("h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts")) { f->Close(); gSystem->Exit(3); }
f->Close(); gSystem->Exit(0);
EOF
[[ $? -eq 0 ]] || fail "RDF crossx output ${rdf_out} has no h2d_crossx_pair_pt_pair_eta_binned_w_signal_cuts — the event loop threw (ROOT exits 0 anyway). Check the log for 'runtime_error'."

# ------ Stage 6: Crossx plotting ------
# use_pt_bins_150=true also refreshes the opt-in pp24_pt_150/ variant. It is NOT
# a second binning of the nominal result (see .claude/CLAUDE.md §Binnings): it is
# the pT_bins_150 axis, plotted so the reach of the data above 120 GeV is visible.
# It must be regenerated in the SAME run as the nominal, or the two directories
# drift apart — they did, from 2026-04-17 to 2026-08-04, because this stage only
# ever called the default (false).
log "Running crossx plotting for PP 2024 (nominal + pt_150 variant)"
pushd "$PLOT_DIR" >/dev/null
root -l -b -q 'plot_single_b_crossx_pp.cxx(24,"",true)'
# RAW pair counts in the same cells as the crossx figure above (PNG + CSV). It reads the
# unweighted twin h2d_counts_pair_pt_pair_eta_binned_w_signal_cuts written by Stage 5, so it must
# run in the same pass -- the statistical reach and the cross-section it belongs to have to come
# from one RDF output. See docs/tracking/pp24_stats_and_powheg_fullsim_compr.md.
root -l -b -q 'plot_pp_counts_pair_pt_in_eta.cxx+(24)'
popd >/dev/null

# ------ Stage 7: Trigger efficiency correction sanity check ------
log "Running before/after trigger efficiency correction sanity plots"
pushd "$PLOT_DIR" >/dev/null
# include_pbpb=false: pp only. The Pb+Pb histograms still carry the retyped 44-bin pair-eta
# axis, on which the coarse panel edge 2.2 is not a bin edge, so PairEtaPanels::Bins throws --
# and it is the FIRST spec, so the pp panels would never be drawn. See the macro's own comment
# and docs/tracking/pp24_all_vertex_pairs.md D5. Restore the default once Pb+Pb is migrated.
root -l -b -q 'plot_crossx_trig_corr_sanity.C(false)'
popd >/dev/null

# ------ Stage 8 (optional): MC-data comparison plots ------
SKIP_MC_DATA_COMPR="${SKIP_MC_DATA_COMPR:-0}"
if [[ "$SKIP_MC_DATA_COMPR" -eq 1 ]]; then
  log "SKIP_MC_DATA_COMPR=1 — skipping MC-data comparison plots"
else
  MC_COMPR_DIR="${ANALYSIS_DIR}/plotting_codes/mc_data_compr"
  log "Running MC-data comparison plots (Pythia pp24 fullsim / POWHEG vs pp data)"
  pushd "$MC_COMPR_DIR" >/dev/null
  # THREE macros since 2026-08-25 (docs/tracking/mc_data_compr_signal_generic_split.md): the set
  # is split into a SIGNAL family (data-like signal-region cuts, so every correction is applied
  # inside the region it was measured in) and a GENERIC family (no signal-region cuts, so the
  # corrections are extrapolated -- see the README in that directory). The signal family is the
  # one that carries physics; running only the generic macro, as this stage used to, would leave
  # the signal plots stale after a crossx rerun.
  root -l -b -q 'plot_mc_data_compr_signal.cxx+()'    # -> plots/mc_data_compr/signal/
  root -l -b -q 'plot_mc_data_pair_pt_in_eta.cxx+()'  # -> plots/mc_data_compr/signal/ (was manual-only)
  root -l -b -q 'plot_mc_data_compr.cxx+()'           # -> plots/mc_data_compr/generic/
  popd >/dev/null
  log "MC-data comparison plots saved to /usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/mc_data_compr/{signal,generic}/"
fi

log "PP crossx/nominal pipeline completed successfully"
