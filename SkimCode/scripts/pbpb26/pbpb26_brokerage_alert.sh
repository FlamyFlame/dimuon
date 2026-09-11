#!/bin/bash
# Proactive brokerage check.  Prints ONE line per task that is in a state worth acting on
# (pending / throttled / "no candidates"), naming the *reason* JEDI gave, so a blocked task
# is never mistaken for a slow one.  Silent when everything is merely running.
#
# Reasons seen on this skim and what they mean:
#   -below_min_weight  user's queued-job budget exhausted (userQRem=0.000)
#   -badsite           user's queue at that site above the 5%-of-site-running cap
#   -cache             GPU/ARM queue without the x86_64 release -- harmless, always present
#   -status            site in test/offline state -- harmless
#   throttled ... type=transfer   task's WAN read volume over JEDI's 6000 GB pacing limit
D=/usatlas/u/yuhanguo/usatlasdata/dimuon_data
T=$(mktemp -d)
for t in $(grep -vE '^[[:space:]]*#|^[[:space:]]*$' "$D/sep2026_pbpb26_skim.txt" | awk '{print $1}'); do
  j=$(curl -s --max-time 90 "https://bigpanda.cern.ch/task/$t/?json" 2>/dev/null)
  [[ -z "$j" ]] && continue
  st=$(python3 -c "
import sys,json
try: d=json.load(sys.stdin); print((d.get('task',d)).get('status') or '')
except Exception: print('')" <<<"$j")
  case "$st" in pending|throttled|exhausted|broken|aborted) ;; *) continue ;; esac
  url=$(python3 -c "
import sys,json,re
try:
  d=json.load(sys.stdin); e=(d.get('task',d)).get('errordialog') or ''
  m=re.search(r'href=\"([^\"]+)\"', e); print(m.group(1) if m else '')
except Exception: print('')" <<<"$j")
  msg=$(python3 -c "
import sys,json,re
try:
  d=json.load(sys.stdin); e=re.sub('<[^>]+>','',(d.get('task',d)).get('errordialog') or '')
  print(e.strip().lstrip(':').strip()[:160])
except Exception: print('')" <<<"$j")
  reasons=""
  if [[ -n "$url" ]]; then
    if curl -s --max-time 90 -o "$T/$t.gz" "$url" 2>/dev/null && gunzip -c "$T/$t.gz" > "$T/$t.txt" 2>/dev/null; then
      # only report reasons that indicate a real block, not the permanent GPU/ARM/test noise
      reasons=$(grep -oE "criteria=-(below_min_weight|badsite|scratch|storage|walltime|pilot|hospital|cap)" "$T/$t.txt" \
                | sort | uniq -c | awk '{printf "%s x%s; ", $2, $1}')
    fi
  fi
  echo "BROKERAGE $t [$st] ${msg} ${reasons:+| blocking reasons: $reasons}"
done
rm -rf "$T"
