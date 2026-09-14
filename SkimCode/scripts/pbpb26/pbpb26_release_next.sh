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
# Health of the newest task, expressed RELATIVELY.  An absolute job-count floor was wrong:
# a recovery task is legitimately small (part6 = 2367 files -> 67 jobs) and can never reach
# a fixed threshold, so the gate blocked part 5 indefinitely.  What actually matters is
# whether the newest task is past scouting and not sitting on a queue backlog.
NEWEST_MAX_ACT_FRAC=70 # activated must be < this % of the task's jobs
COOLDOWN=5400          # 90 min
MAX_LIVE=3             # release only when at most this many tasks are still live.
                       # Part 5's runs (523023-523437) have DATADISK replicas at only
                       # RAL and BNL; with four other tasks running, the user is over
                       # the 5%-of-site-running queue cap at BOTH, so part 5 could not
                       # be brokered anywhere and sat at 37/20334 files.  Waiting for a
                       # task to finish is the only thing that frees those two sites.

log() { echo "[$(date -u +%FT%TZ)] $*" >> "$LOG"; }

next=$(grep -vE '^\s*#|^\s*$' "$PENDING" 2>/dev/null | head -1)
[[ -z "$next" ]] && exit 0                     # nothing left to release

# Independent guard against double-submission: if a task for this part already exists in
# the bookkeeping file, drop it from the pending list and stop.  This catches any future
# bookkeeping slip without relying on the pending list alone.
if grep -qE "PbPb2026data\.[^ ]*\.part${next}\._EXT0" "$BOOK" 2>/dev/null; then
  log "part $next already submitted (present in $BOOK) -- removing from pending list"
  grep -vE "^[[:space:]]*${next}[[:space:]]*$" "$PENDING" > "$PENDING.tmp" || true
  mv "$PENDING.tmp" "$PENDING"
  exit 0
fi

now=$(date +%s)
last=$(cat "$STAMP" 2>/dev/null || echo 0)
(( now - last < COOLDOWN )) && exit 0

tasks=$(grep -vE '^\s*#|^\s*$' "$BOOK" | awk '{print $1}')
nlive=0
for t in $tasks; do
  tst=$(curl -s --max-time 60 "https://bigpanda.cern.ch/task/$t/?json" 2>/dev/null | python3 -c "
import sys,json
try: d=json.load(sys.stdin); print((d.get('task',d)).get('status') or '')
except Exception: print('')" 2>/dev/null)
  case "$tst" in done|finished|failed|broken|aborted) ;; *) nlive=$((nlive+1)) ;; esac
done
if (( nlive > MAX_LIVE )); then
  log "hold part $next: $nlive tasks still live (> $MAX_LIVE)"
  exit 0
fi
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
newest_st=$(curl -s --max-time 60 "https://bigpanda.cern.ch/task/$newest/?json" 2>/dev/null | python3 -c "
import sys,json
try: d=json.load(sys.stdin); print((d.get('task',d)).get('status') or '')
except Exception: print('')" 2>/dev/null)
case "$newest_st" in
  scouting|pending|registered|defined) ;;   # not yet proven to broker -- hold
  *)
    if (( newest_jobs > 0 && newest_prog > 0 \
          && newest_act * 100 < newest_jobs * NEWEST_MAX_ACT_FRAC )); then
      newest_ok=1
    fi ;;
esac

if (( total_act >= ACT_MAX )); then
  log "hold part $next: total activated=$total_act >= $ACT_MAX"
  exit 0
fi
if (( newest_ok == 0 )); then
  log "hold part $next: newest task $newest not scaled/running yet (jobs=$newest_jobs prog=$newest_prog act=$newest_act; status=$newest_st; need non-scouting, prog>0, act<${NEWEST_MAX_ACT_FRAC}% of jobs)"
  exit 0
fi

log "RELEASING part $next (total activated=$total_act, newest task $newest healthy)"
cd /usatlas/u/yuhanguo/workarea/dimuon_codes/SkimCode || exit 1
set +u
source setup_26.sh  >> "$LOG" 2>&1
lsetup panda        >> "$LOG" 2>&1
set -u
cd run_26hi || exit 1
# Prefer a _v3 variant when one exists.  pathena REACTIVATES an existing task if the outDS
# tag is unchanged -- "reactivation accepted ... will be re-executed" -- and a reactivated
# task keeps its ORIGINAL sandbox, i.e. the pre-ZDC-fix library.  A bumped tag forces a new
# task with the current build.
script="grid_sub_part${next}.sh"
[[ -f "grid_sub_part${next}_v3.sh" ]] && script="grid_sub_part${next}_v3.sh"
out=$(bash "$script" 2>&1); echo "$out" >> "$LOG"
tid=$(grep -oE 'new jediTaskID=[0-9]+' <<<"$out" | grep -oE '[0-9]+' | head -1)
if [[ -z "$tid" ]]; then
  log "SUBMISSION FAILED for part $next -- left in the pending list"
  echo "PbPb2026: submission of part $next FAILED (see $LOG)"
  exit 1
fi
echo "$tid user.yuhang.TrigRates.dimuon.PbPb2026data.Sep2026.$(grep -oE "Sep2026\\.v[0-9]+" "$script" | head -1 | cut -d. -f2).part${next}._EXT0" >> "$BOOK"
# NOTE: `grep -v ... > tmp && mv` is WRONG here.  When the removed entry was the LAST
# one, grep prints nothing and exits 1, so the `&&` short-circuits and the pending list is
# never truncated -- the part stays queued and gets released a SECOND time once the
# cooldown expires, producing a duplicate task and a second output dataset that would be
# hadded into the same data_pbpb26_part<N>.root.  Ignore grep's exit status.
grep -vE "^[[:space:]]*${next}[[:space:]]*$" "$PENDING" > "$PENDING.tmp" || true
mv "$PENDING.tmp" "$PENDING"
date +%s > "$STAMP"
log "part $next submitted as jediTaskID=$tid"
echo "PbPb2026: released part $next -> jediTaskID=$tid (total activated was $total_act)"
