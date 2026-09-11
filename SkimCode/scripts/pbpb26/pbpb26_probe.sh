#!/bin/bash
# One probe of the 5 PbPb2026 skim tasks.  Prints one summary line per task at the JOB
# level (task-level nfilesfinished lags badly: it only counts files whose MERGE finished,
# so a task can look frozen for hours while thousands of jobs are merging).
# Appends to the progress log and echoes any line worth acting on.
BOOK=/usatlas/u/yuhanguo/usatlasdata/dimuon_data/sep2026_pbpb26_skim.txt
# Read the live task IDs from the bookkeeping file so a newly released part is picked up
# automatically, with no edit here.
TASKS=$(grep -vE "^\s*#|^\s*$" "$BOOK" | awk '{print $1}' | tr '\n' ' ')
LOG=/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb26_task_progress.log
stamp=$(date -u +%FT%TZ)
for t in $TASKS; do
  curl -s --max-time 90 "https://bigpanda.cern.ch/jobs/?jeditaskid=$t&json&limit=20000" 2>/dev/null | python3 -c "
import sys,json,collections
try:
    d=json.load(sys.stdin); jobs=d.get('jobs',[])
    c=collections.Counter(j.get('jobstatus') for j in jobs)
    done=c.get('finished',0); merg=c.get('merging',0); fail=c.get('failed',0)
    act=c.get('activated',0); run=c.get('running',0)+c.get('starting',0)
    print(f\"$t jobs={len(jobs)} fin={done} merg={merg} run={run} act={act} fail={fail}\")
except Exception:
    print(f\"$t QUERYFAIL\")
" 2>/dev/null
done | tee -a "$LOG" | sed "s/^/$stamp /"
