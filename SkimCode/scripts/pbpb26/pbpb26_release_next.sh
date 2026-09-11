#!/bin/bash
# Release the next held PbPb2026 part when the PanDA queued-job budget has freed up.
#
# The v1 submission failed because all five parts were released at once: ~8300 queued jobs
# drove the user's brokerage weight below MIN_WEIGHT_user (userQRem=0.000), which blocked
# the MERGE brokerage of every task.  So parts are released one at a time, and only when:
#   (1) the most recently released task has SCALED UP and is demonstrably running --
#       >= NEWEST_MIN_JOBS jobs, >= NEWEST_MIN_PROG of them running/merging/finished, and
#       activated not the overwhelming majority.  A just-submitted task holds only a
#       handful of scout jobs, so "it has some jobs" is NOT evidence of spare capacity;
#       waiting for it to scale is.
#   (2) total `activated` across live tasks below ACT_MAX -- a BACKSTOP only.  A raw
#       activated count is a poor jam signal by itself: the v1 jam happened at ~8300
#       queued, yet 52488080 later ran happily at 8162 queued once EMMY_KIT stopped
#       hoarding.  What actually jammed was ~4500 jobs parked at a site running none of
#       them -- which v2 avoids by excluding that site.
#   (3) at least COOLDOWN seconds since the last release.
# Prints one line when it releases (or when it declines for a reason worth seeing).
set -uo pipefail
D=/usatlas/u/yuhanguo/usatlasdata/dimuon_data
BOOK=$D/sep2026_pbpb26_skim.txt
PENDING=$D/pbpb26_pending_parts.txt
STAMP=$D/pbpb26_last_release.stamp
LOG=$D/pbpb26_release.log
ACT_MAX=15000          # backstop only -- see the note above
NEWEST_MIN_JOBS=200    # the newest task must have scaled past its scout phase
NEWEST_MIN_PROG=100    # ...and actually be executing, not merely queued
COOLDOWN=5400          # 90 min

log() { echo "[$(date -u +%FT%TZ)] $*" >> "$LOG"; }

next=$(grep -vE '^\s*#|^\s*$' "$PENDING" 2>/dev/null | head -1)
[[ -z "$next" ]] && exit 0                     # nothing left to release

now=$(date +%s)
last=$(cat "$STAMP" 2>/dev/null || echo 0)
(( now - last < COOLDOWN )) && exit 0

tasks=$(grep -vE '^\s*#|^\s*$' "$BOOK" | awk '{print $1}')
# The "newest" task is the one with the LARGEST jediTaskID -- PanDA hands them out
# monotonically.  Do NOT use the last line of the bookkeeping file: the v1 survivor
# 52488080 sits below the v2 part-1 task 52491225 there, so "last line" picked an OLDER,
# already-scaled task and the scale-up gate passed vacuously (this released part 3 early).
total_act=0; newest=""; newest_jobs=0; newest_prog=0; newest_act=0
for t in $tasks; do
  r=$(curl -s --max-time 90 "https://bigpanda.cern.ch/jobs/?jeditaskid=$t&json&limit=30000" 2>/dev/null | python3 -c "
import sys,json,collections
try:
    d=json.load(sys.stdin); j=d.get('jobs',[])
    c=collections.Counter(x.get('jobstatus') for x in j)
    print(c.get('activated',0), len(j), c.get('running',0)+c.get('merging',0)+c.get('finished',0)+c.get('starting',0))
except Exception: print(-1,-1,-1)
" 2>/dev/null)
  read -r act njobs prog <<<"$r"
  [[ "$act" == "-1" || -z "$act" ]] && { log "probe failed for $t -- not releasing this cycle"; exit 0; }
  total_act=$(( total_act + act ))
  if [[ -z "$newest" || "$t" -gt "$newest" ]]; then
    newest="$t"; newest_jobs=$njobs; newest_prog=$prog; newest_act=$act
  fi
done
newest_ok=0
if (( newest_jobs >= NEWEST_MIN_JOBS && newest_prog >= NEWEST_MIN_PROG \
      && newest_act * 10 < newest_jobs * 7 )); then
  newest_ok=1
fi

if (( total_act >= ACT_MAX )); then
  log "hold part $next: total activated=$total_act >= $ACT_MAX"
  exit 0
fi
if (( newest_ok == 0 )); then
  log "hold part $next: newest task $newest not scaled/running yet (jobs=$newest_jobs prog=$newest_prog act=$newest_act; need jobs>=$NEWEST_MIN_JOBS prog>=$NEWEST_MIN_PROG act<70%)"
  exit 0
fi

log "RELEASING part $next (total activated=$total_act, newest task $newest healthy)"
cd /usatlas/u/yuhanguo/workarea/dimuon_codes/SkimCode || exit 1
set +u
source setup_26.sh  >> "$LOG" 2>&1
lsetup panda        >> "$LOG" 2>&1
set -u
cd run_26hi || exit 1
out=$(bash "grid_sub_part${next}.sh" 2>&1); echo "$out" >> "$LOG"
tid=$(grep -oE 'new jediTaskID=[0-9]+' <<<"$out" | grep -oE '[0-9]+' | head -1)
if [[ -z "$tid" ]]; then
  log "SUBMISSION FAILED for part $next -- left in the pending list"
  echo "PbPb2026: submission of part $next FAILED (see $LOG)"
  exit 1
fi
echo "$tid user.yuhang.TrigRates.dimuon.PbPb2026data.Sep2026.v2.part${next}._EXT0" >> "$BOOK"
grep -vE "^\s*${next}\s*$" "$PENDING" > "$PENDING.tmp" && mv "$PENDING.tmp" "$PENDING"
date +%s > "$STAMP"
log "part $next submitted as jediTaskID=$tid"
echo "PbPb2026: released part $next -> jediTaskID=$tid (total activated was $total_act)"
