#!/bin/bash
# Print one line per PbPb2026 task, flagging the WEDGED signature:
#   live jobs (running+activated+defined+starting+assigned) == 0
#   AND a non-empty merging backlog
#   AND that backlog untouched for > COLD_H hours
# This is the only pattern seen on this campaign that a user can neither wait out nor fix
# with retryTask.  Everything else -- throttled pacing, bursty terminal counts, tasks
# flipping pending/running -- looks like a stall but is not one.
COLD_H=${COLD_H:-6}
D=/usatlas/u/yuhanguo/usatlasdata/dimuon_data
for t in $(grep -vE '^[[:space:]]*#|^[[:space:]]*$' "$D/sep2026_pbpb26_skim.txt" | awk '{print $1}'); do
  curl -s --max-time 90 "https://bigpanda.cern.ch/jobs/?jeditaskid=$t&json&limit=40000" 2>/dev/null | \
  COLD_H="$COLD_H" TID="$t" python3 -c "
import sys,json,collections,datetime,os
try:
    d=json.load(sys.stdin); j=d.get('jobs',[])
except Exception:
    print(f\"{os.environ['TID']} QUERYFAIL\"); raise SystemExit
c=collections.Counter(x.get('jobstatus') for x in j)
live=sum(c.get(k,0) for k in ('running','activated','defined','starting','assigned'))
mg=[x for x in j if x.get('jobstatus')=='merging']
age=None
if mg:
    newest=max((x.get('modificationtime') or '') for x in mg)
    try: age=(datetime.datetime.utcnow()-datetime.datetime.fromisoformat(newest)).total_seconds()/3600
    except Exception: pass
cold=float(os.environ['COLD_H'])
wedged = (live==0 and mg and age is not None and age>cold)
print(f\"{os.environ['TID']} live={live} merging={len(mg)} merge_cold={'%.1fh'%age if age is not None else '-'}\" + (' WEDGED' if wedged else ''))
"
done
