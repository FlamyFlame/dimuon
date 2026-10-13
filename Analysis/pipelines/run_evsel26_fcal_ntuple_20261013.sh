#!/usr/bin/env bash
# One-off (2026-10-13), tracking doc docs/tracking/ntuple_rerun_pp24_pbpb2526_20261013.md:
#  A. re-derive the Pb+Pb 2026 event selection with D8 (522200/522041 dropped) + cut-2 window, in the standard
#     order: event_level (cut 1, nominal + alt band) -> cuts (2-5) -> cuts_alt (2-5);
#  B. wait for the 2025 cut-2 chain of another session (PID $P25CHAIN) to end, require fresh 2025 cuts files;
#  C. the 4 FCal comparisons (compiled once, then run in parallel, distinct outputs);
#  D. Condor NTuple reruns pbpb25 + pbpb26 (run_ntuple_rerun_20261013.sh).
# Every ROOT step runs under a wall-clock deadline WITHOUT blocking on the process (a /pnfs NFS read in D state
# cannot be killed and would block `wait`/`timeout`): on deadline -> FATAL and exit. Last line: ALL DONE / FATAL.
source /usatlas/u/yuhanguo/setup.sh >/dev/null 2>&1
export X509_USER_PROXY=/tmp/x509up_u101379
set -uo pipefail
A=/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis
P25CHAIN=${P25CHAIN:?pid of the other session run_cut2_chain.sh}
T25=${T25:?epoch the 2025 chain started}
L=$A/pipelines/logs/evsel26_fcal_20261013; mkdir -p $L
BD=$L/build; mkdir -p $BD    # own ACLiC build dir: never touches the other session's .so
T0=$(date +%s)
log(){ echo "[$(date '+%F %T')] $*"; }
fail(){ log "FATAL: $*"; exit 1; }

# compile_once <macro.cxx>: ACLiC into $BD, error code checked
compile_once(){
  root -l -b -q -e "gSystem->SetBuildDir(\"$BD\", true); int e=0; gROOT->ProcessLine(\".L $1+\", &e); gSystem->Exit(e ? 2 : 0);" \
    > $L/compile_$(basename $1 .cxx).log 2>&1 < /dev/null || fail "compile $1 (see $L)"
}
# run_dl <deadline_s> <name> <macro.cxx> <call>: load the compiled macro, run the call; any error code -> rc 2/3
run_dl(){
  local dl=$1 name=$2 mac=$3 call=$4 t
  t=$(date +%s)
  cat > $L/$name.C <<EOF
{
  gSystem->SetBuildDir("$BD", true);
  int e = 0;
  gROOT->ProcessLine(".L $mac+", &e);
  if (e) gSystem->Exit(2);
  gROOT->ProcessLine("$call", &e);
  gSystem->Exit(e ? 3 : 0);
}
EOF
  root -l -b -q $L/$name.C > $L/$name.log 2>&1 < /dev/null &
  local pid=$!
  while kill -0 $pid 2>/dev/null; do
    if (( $(date +%s) - t > dl )); then kill $pid 2>/dev/null; fail "$name exceeded ${dl}s (hang?) - see $L/$name.log"; fi
    sleep 20
  done
  wait $pid || fail "$name rc=$? - see $L/$name.log"
  grep -qE '^Error in|terminate called|Segmentation' $L/$name.log && fail "$name: error lines in $L/$name.log"
  log "$name OK ($(( ($(date +%s)-t)/60 )) min)"
}

D26=~/usatlasdata/dimuon_data/pbpb_2026; D25=~/usatlasdata/dimuon_data/pbpb_2025
cd $A/plotting_codes/event_selection
log "A: event_level(26) -> cuts(26) -> cuts_alt(26)"
compile_once plot_pbpb_event_sel_event_level.cxx
compile_once plot_pbpb_event_sel_cuts.cxx
run_dl 7200 event_level_26 plot_pbpb_event_sel_event_level.cxx "plot_pbpb_event_sel_event_level(26)"
run_dl 7200 cuts_26        plot_pbpb_event_sel_cuts.cxx        "plot_pbpb_event_sel_cuts(26)"
run_dl 7200 cuts_alt_26    plot_pbpb_event_sel_cuts.cxx        "plot_pbpb_event_sel_cuts_alt(26)"
for f in $D26/event_sel_cuts_pbpb_2026.root $D26/event_sel_cuts_pbpb_2026_alt.root; do
  (( $(stat -c %Y $f) > T0 )) || fail "$f not refreshed"
done
for n in event_level_26 cuts_26 cuts_alt_26; do
  grep -q 'Excluding 4 bad run(s)' $L/$n.log || fail "$n did not apply the 4 D7+D8 bad runs"
done

log "B: waiting for the 2025 chain (pid $P25CHAIN)"
tw=$(date +%s)
while kill -0 $P25CHAIN 2>/dev/null; do
  (( $(date +%s) - tw > 10800 )) && fail "2025 chain still running after 3 h"
  sleep 60
done
for f in $D25/event_sel_cuts_pbpb_2025.root $D25/event_sel_cuts_pbpb_2025_alt.root; do
  (( $(stat -c %Y $f) > T25 )) || fail "$f not refreshed by the 2025 chain"
done
log "2025 cuts files fresh"

log "C: 4 FCal comparisons in parallel"
compile_once plot_pbpb_fcal_comparison.cxx
pids=()
for fn in plot_pbpb_fcal_comparison plot_pbpb_fcal_comparison_alt plot_pbpb_fcal_comparison_262524 plot_pbpb_fcal_comparison_262524_alt; do
  ( run_dl 14400 $fn plot_pbpb_fcal_comparison.cxx "$fn()"
    grep -q 'png file .* has been created' $L/$fn.log || fail "$fn wrote no PNG" ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait $p || fail "an FCal comparison failed (see $L)"; done

log "D: NTuple reruns pbpb25 + pbpb26"
cd $A/pipelines
CUTS_AFTER=$T25 ./run_ntuple_rerun_20261013.sh pbpb25 > logs/ntuple_rerun_20261013_pbpb25.log 2>&1 & q25=$!
CUTS_AFTER=$T0  ./run_ntuple_rerun_20261013.sh pbpb26 > logs/ntuple_rerun_20261013_pbpb26.log 2>&1 & q26=$!
wait $q25; r25=$?; wait $q26; r26=$?
tail -1 logs/ntuple_rerun_20261013_pbpb25.log; tail -1 logs/ntuple_rerun_20261013_pbpb26.log
(( r25 == 0 && r26 == 0 )) || fail "NTuple rerun rc pbpb25=$r25 pbpb26=$r26"
log "ALL DONE"
