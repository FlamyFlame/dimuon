#!/bin/bash
# For every live PbPb2026 task: report status + errordialog, and whenever the task is
# pending / brokerage-failing, fetch and decode its JEDI brokerage log and print the
# candidate-cut summary and the reason each surviving site was rejected.
#
# The jedilog host is per-task and is embedded in the task's errordialog -- do not
# hardcode aipanda094/095.  The file is gzipped.
D=/usatlas/u/yuhanguo/usatlasdata/dimuon_data
T=${TMPDIR:-/tmp}/pbpb26_brok.$$
mkdir -p "$T"
for t in $(grep -vE '^[[:space:]]*#|^[[:space:]]*$' "$D/sep2026_pbpb26_skim.txt" | awk '{print $1}'); do
  j=$(curl -s --max-time 90 "https://bigpanda.cern.ch/task/$t/?json" 2>/dev/null)
  read -r st err <<<"$(python3 -c "
import sys,json,re
try:
    d=json.load(sys.stdin); k=d.get('task',d)
    e=(k.get('errordialog') or '-')
    url=''
    m=re.search(r'href=\"([^\"]+)\"', e)
    if m: url=m.group(1)
    e=re.sub('<[^>]+>','',e).strip().lstrip(':').strip()
    print(k.get('status'), url or '-')
except Exception: print('?','-')
" <<<"$j")"
  echo "=============== task $t : status=$st"
  python3 -c "
import sys,json,re
d=json.load(sys.stdin); k=d.get('task',d)
e=re.sub('<[^>]+>','',(k.get('errordialog') or '-')).strip().lstrip(':').strip()
print('   errordialog:', e if e else '(none)')
" <<<"$j"
  [[ "$err" == "-" ]] && continue
  curl -s --max-time 120 -o "$T/$t.gz" "$err" 2>/dev/null || continue
  gunzip -c "$T/$t.gz" > "$T/$t.txt" 2>/dev/null || { echo "   (jedilog not decodable/expired)"; continue; }
  echo "   --- brokerage summary ---"
  sed -n '/Job brokerage summary/,/no candidates\|final candidates/p' "$T/$t.txt" | sed 's/^[0-9-]* [0-9:.]* : /     /'
  echo "   --- distinct rejection reasons ---"
  grep -oE "criteria=-[a-z_]+" "$T/$t.txt" | sort | uniq -c | sed 's/^/     /'
  echo "   --- one example of each reason ---"
  for c in $(grep -oE "criteria=-[a-z_]+" "$T/$t.txt" | sort -u); do
    grep -m1 -- "$c" "$T/$t.txt" | sed 's/^[0-9-]* [0-9:.]* : /     /' | cut -c1-300
  done
done
rm -rf "$T"
