#!/usr/bin/env bash
# =============================================================================
# Round-7: re-run the store_mc_trigger SINGLE-MUON NTP passes so the trees carry the new
# per-muon `n_vtx` (number of reconstructed, track-bearing primary vertices in the event;
# MuonFullsimExtra::n_vtx, filled from the skim's vtx_ntrk). Required by the Step-1
# SANITY CHECK (mc_trigger_efficiency.md §3.5, contract item 2).
#
# ONLY the single-muon trees are redone -- the sanity check is a Step-1 study, and nothing
# else in the round-7 product set needs a vertex count. The pair trees are NOT re-run.
#
# Everything else in these passes is unchanged, so the trees must come back with IDENTICAL
# entry counts; the script checks that against the pre-run values.
#
# Usage: ./run_ntp_singles_round7_nvtx.sh          (all three samples)
#        SAMPLES="overlay noovl" ./run_ntp_singles_round7_nvtx.sh
# =============================================================================
set -Euo pipefail

SAMPLES="${SAMPLES:-pp_full overlay noovl}"
ANALYSIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NTP_DIR="${ANALYSIS_DIR}/NTupleProcessingCode"
LOG_DIR="${ANALYSIS_DIR}/pipelines/logs_round7"
mkdir -p "$LOG_DIR"

ts(){ date "+%Y-%m-%d %H:%M:%S"; }
log(){ echo "[$(ts)] $*"; }
fail(){ echo "[$(ts)] FATAL: $*" >&2; exit 1; }

script_of(){ case "$1" in
  pp_full) echo run_pythia_fullsim_single_muon_mc_trig_full_sample.sh ;;
  overlay) echo run_pythia_fullsim_overlay_single_muon_mc_trig.sh ;;
  noovl)   echo run_pythia_fullsim_noovl_single_muon_mc_trig.sh ;;
  *) fail "unknown sample $1" ;; esac; }
# Sample dirs from the shell twin of FullSimSampleType.h. The round-7 n_vtx trees were made on
# the pbpb23 overlay (r17618), so this script pins OVERLAY_YEAR=23 unless told otherwise.
: "${OVERLAY_YEAR:=23}"
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/fullsim_sample_layout.sh"
tree_of(){ case "$1" in
  pp_full|overlay|noovl)
    echo "$(fullsim_sample_dir "$1")muon_pairs_pythia_fullsim_$(fullsim_sample_label "$1")_no_data_resonance_cuts_mc_trig_single_muon$([[ $1 == pp_full ]] && echo _full).root" ;;
  *) fail "unknown sample $1" ;; esac; }

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
set +eu
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh" --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"
set -eu

entries(){ root -l -b -q -e "TFile*f=TFile::Open(\"$1\"); TTree*t=f&&!f->IsZombie()?(TTree*)f->Get(\"muon_tree\"):nullptr; printf(\"N %lld\n\", t?t->GetEntries():-1LL);" 2>/dev/null | grep -oP 'N \K-?[0-9]+' | tail -1; }
has_nvtx(){ root -l -b -q -e "TFile*f=TFile::Open(\"$1\"); TTree*t=(TTree*)f->Get(\"muon_tree\"); printf(\"HAS %d\n\", (t&&t->GetLeaf(\"n_vtx\"))?1:0);" 2>/dev/null | grep -oP 'HAS \K[0-9]+' | tail -1; }

# record the pre-run entry counts (the re-run must reproduce them exactly)
declare -A BEFORE
for s in $SAMPLES; do BEFORE[$s]=$(entries "$(tree_of "$s")"); log "before: $s = ${BEFORE[$s]} entries"; done

# pre-compile ONCE (concurrent ACLiC builds race on the shared _h.so)
log "pre-compiling PythiaAnalysisClasses.h+"
( cd "$NTP_DIR" && root -l -b -q -e '.L PythiaAnalysisClasses.h+' ) >"$LOG_DIR/compile_ntp.log" 2>&1 \
    || fail "PythiaAnalysisClasses compile failed (see $LOG_DIR/compile_ntp.log)"

pids=(); names=()
for s in $SAMPLES; do
    log "launching NTP singles: $s -> $LOG_DIR/ntp_${s}.log"
    ( cd "$NTP_DIR" && bash "$(script_of "$s")" ) >"$LOG_DIR/ntp_${s}.log" 2>&1 &
    pids+=($!); names+=("$s")
done
rc=0
for i in "${!pids[@]}"; do
    if wait "${pids[$i]}"; then log "  NTP OK : ${names[$i]}"; else log "  NTP FAILED : ${names[$i]}"; rc=1; fi
done
[[ $rc -eq 0 ]] || fail "an NTP pass failed -- see $LOG_DIR/ntp_*.log"

# a ROOT macro that throws still exits 0 -> validate the artefacts
for s in $SAMPLES; do
    f="$(tree_of "$s")"
    n=$(entries "$f"); hv=$(has_nvtx "$f")
    [[ "$hv" == "1" ]] || fail "$s: n_vtx leaf still missing in $f"
    [[ "$n" == "${BEFORE[$s]}" ]] || fail "$s: entry count changed ${BEFORE[$s]} -> $n (the vertex propagation must not change the sample)"
    log "  validated $s : $n entries, n_vtx present"
done
log "════ round-7 singles NTP (n_vtx) DONE ════"
