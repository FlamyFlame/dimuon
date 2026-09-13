#!/bin/bash
# Force grid_monitor to RE-download a task whose input-file count has grown since it was
# merged.
#
# Why this is needed: a task can reach terminal state 'finished' (not 'done') with input
# files left unprocessed, and grid_monitor happily downloads and merges that partial output
# -- it accepts "done or finished (>=90% file success)".  If the task is then retried and
# recovers those files, grid_monitor has already marked it 'completed' and will never look
# again, so the extra events would be silently missing from data_pbpb26_part<N>.root.
#
# Strategy: snapshot nfilesfinished when a task is first seen as 'completed'.  If it later
# exceeds the snapshot, reset the state to 'pending' and clear the snapshot -- grid_monitor
# re-downloads (rucio skips files already on disk) and re-hadds.
D=/usatlas/u/yuhanguo/usatlasdata/dimuon_data
SNAP=$D/pbpb26_downloaded_snapshot.txt
STATE=$D/grid_monitor_state.txt
touch "$SNAP"
for t in $(grep -vE '^[[:space:]]*#|^[[:space:]]*$' "$D/sep2026_pbpb26_skim.txt" | awk '{print $1}'); do
  st=$(grep "^${t} " "$STATE" 2>/dev/null | head -1 | awk '{print $2}')
  fin=$(curl -s --max-time 90 "https://bigpanda.cern.ch/task/$t/?json" 2>/dev/null | python3 -c "
import sys,json
try:
    d=json.load(sys.stdin); ds=d.get('datasets',[])
    print(sum(x.get('nfilesfinished',0) for x in ds if x.get('type')=='input'))
except Exception: print(-1)" 2>/dev/null)
  [[ "$fin" == "-1" || -z "$fin" ]] && continue
  if [[ "$st" == "completed" ]]; then
    prev=$(grep "^${t} " "$SNAP" | head -1 | awk '{print $2}')
    if [[ -z "$prev" ]]; then
      echo "${t} ${fin}" >> "$SNAP"
    elif (( fin > prev )); then
      # Only re-download once the recovery has SETTLED.  A retried task trickles files back
      # in over hours; resetting on the first increment makes grid_monitor re-download and
      # re-merge ~50 GB for a merge that is still partial, over and over.  Require the task
      # to be in a terminal state AND its file count unchanged since the previous poll.
      st_panda=$(curl -s --max-time 90 "https://bigpanda.cern.ch/task/$t/?json" 2>/dev/null | python3 -c "
import sys,json
try: d=json.load(sys.stdin); print((d.get('task',d)).get('status') or '')
except Exception: print('')" 2>/dev/null)
      last=$(grep "^${t} " "$SNAP.pending" 2>/dev/null | head -1 | awk '{print $2}')
      sed -i "/^${t} /d" "$SNAP.pending" 2>/dev/null; echo "${t} ${fin}" >> "$SNAP.pending"
      case "$st_panda" in
        done|finished)
          if [[ "$last" == "$fin" ]]; then
            sed -i "s/^${t} completed.*/${t} pending/" "$STATE"
            sed -i "/^${t} /d" "$SNAP"; sed -i "/^${t} /d" "$SNAP.pending" 2>/dev/null
            echo "RE-DOWNLOAD queued: task $t recovered $((fin-prev)) input files and has settled at $fin; grid_monitor state reset to pending"
          fi ;;
      esac
    fi
  fi
done
