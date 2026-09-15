#!/bin/bash
# Standalone PbPb2026 watcher.  Runs under setsid so it survives the Claude session that
# started it -- a session restart on 2026-09-15 killed every in-session monitor, and with
# them the one that was supposed to restart grid_monitor after the dCache outage.  All
# output goes to a log the session can tail; nothing here depends on the session living.
D=/usatlas/u/yuhanguo/usatlasdata/dimuon_data
LOG=$D/pbpb26_watcher.log
YEAR_FILES=111279
prevkey=""
declare -A wedged incomplete
say() { echo "[$(date -u +%FT%TZ)] $*" >> "$LOG"; }
say "watcher started pid=$$"
while true; do
  gp=$(cat $D/pbpb26_grid_monitor.pid 2>/dev/null)
  if [ -z "$gp" ] || ! kill -0 "$gp" 2>/dev/null; then
    # Only restart while there is something left to download.  grid_monitor exits by
    # itself once every task in its list is terminal ("All tasks resolved"); restarting it
    # then just makes it exit again every cycle.  A task not yet 'completed'/'failed' in
    # the state file, or a task in the bookkeeping file with no state entry (newly
    # submitted), means there IS work.
    if grep -vE '^[[:space:]]*#|^[[:space:]]*$' $D/sep2026_pbpb26_skim.txt | awk '{print $1}' | while read -r t; do
         st=$(grep "^$t " $D/grid_monitor_state.txt | awk '{print $2}')
         [ "$st" != "completed" ] && [ "$st" != "failed" ] && echo work
       done | grep -q work; then
      say "WATCHDOG: grid_monitor was dead with work pending -- $(bash $D/pbpb26_grid_monitor_start.sh 2>&1 | tail -1)"
    fi
  fi
  # Second grid_monitor instance for the May-2026 recovery parts (own bookkeeping + pid file).
  mp=$(cat $D/may26rec_grid_monitor.pid 2>/dev/null)
  if [ -z "$mp" ] || ! kill -0 "$mp" 2>/dev/null; then
    if grep -vE '^[[:space:]]*#|^[[:space:]]*$' $D/sep2026_may26_recovery.txt 2>/dev/null | awk '{print $1}' | while read -r t; do
         st=$(grep "^$t " $D/grid_monitor_state.txt | awk '{print $2}')
         [ "$st" != "completed" ] && [ "$st" != "failed" ] && echo work
       done | grep -q work; then
      say "WATCHDOG: may26rec grid_monitor was dead with work pending -- $(bash $D/may26rec_grid_monitor_start.sh 2>&1 | tail -1)"
    fi
  fi
  bash $D/pbpb26_recheck_completed.sh 2>/dev/null | while read -r l; do say "$l"; done
  bash $D/pbpb26_release_next.sh 2>&1 | grep -E "released part|FAILED" | while read -r l; do say "$l"; done
  w=$(bash $D/pbpb26_wedge_check.sh 2>/dev/null)
  while read -r tid rest; do
    [ -z "$tid" ] && continue
    if grep -q "WEDGED" <<<"$rest"; then
      [ "${wedged[$tid]:-0}" = "0" ] && { say "WEDGED: task $tid -- $rest"; wedged[$tid]=1; }
    else
      [ "${wedged[$tid]:-0}" = "1" ] && say "RECOVERED: task $tid -- $rest"; wedged[$tid]=0
    fi
  done <<<"$w"
  snap=$(for t in $(cat $D/sep2026_pbpb26_skim.txt $D/sep2026_may26_recovery.txt 2>/dev/null | grep -vE '^[[:space:]]*#|^[[:space:]]*$' | awk '{print $1}'); do
    curl -s --max-time 60 "https://bigpanda.cern.ch/task/$t/?json" 2>/dev/null | python3 -c "
import sys,json
try:
    d=json.load(sys.stdin); k=d.get('task',d); ds=d.get('datasets',[])
    tot=sum(x.get('nfiles',0) for x in ds if x.get('type')=='input')
    fin=sum(x.get('nfilesfinished',0) for x in ds if x.get('type')=='input')
    print(f\"$t|{k.get('status')}|{fin}|{tot}\")
except Exception: pass" 2>/dev/null
  done)
  while IFS='|' read -r t st fin tot; do
    [ -z "$t" ] && continue
    grep -q "^$t " $D/pbpb26_covered_gaps.txt 2>/dev/null && continue
    if [ "$st" = "finished" ] && [ "${fin:-0}" -lt "${tot:-0}" ]; then
      miss=$((tot-fin))
      [ "${incomplete[$t]:-}" != "$miss" ] && { say "INCOMPLETE: task $t finished $fin/$tot ($miss MISSING, not yet covered)"; incomplete[$t]=$miss; }
    elif [ "$st" = "done" ]; then
      [ "${incomplete[$t]:-}" != "done" ] && { say "TERMINAL: task $t -> done ($fin/$tot, COMPLETE)"; incomplete[$t]=done; }
    fi
  done <<<"$snap"
  key=$(awk -F'|' 'NF>3{printf "%s:%d ", $1, ($4?10*int(10*$3/$4+0.5):0)}' <<<"$snap")
  if [ "$key" != "$prevkey" ]; then
    say "--- progress ---"
    awk -F'|' 'NF>3{printf "%s %-10s %d%%  (%d/%d files)\n", $1, $2, ($4?10*int(10*$3/$4+0.5):0), $3, $4}' <<<"$snap" | while read -r l; do say "  $l"; done
    awk -F'|' -v Y="$YEAR_FILES" 'NF>3{f+=$3} END {printf "YEAR %d/%d = %.1f%%\n", f, Y, 100*f/Y}' <<<"$snap" | while read -r l; do say "  $l"; done
    prevkey="$key"
  fi
  sleep 1800
done
