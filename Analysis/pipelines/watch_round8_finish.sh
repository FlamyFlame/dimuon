#!/usr/bin/env bash
# Waits for the round-8 both-binnings driver to exit, then writes ONE summary file.
# Local on purpose: the state it inspects (GPFS ROOT files, driver PID, stage logs) does not
# exist anywhere but this cluster, so a cloud routine could not see any of it.
A=/gpfs/mnt/atlasgpfs01/usatlas/workarea/yuhanguo/dimuon_codes/Analysis
L=$A/pipelines/logs_round8
OUT=$L/ROUND8_SUMMARY.txt
while pgrep -f run_round8_mc_both_binnings >/dev/null 2>&1; do sleep 60; done
{
  echo "round-8 both-binnings driver finished at $(date '+%F %H:%M')"
  echo
  echo "== driver =="; tail -6 "$L/both_binnings_driver.log" 2>/dev/null
  for m in 8bin 4bin; do
    echo; echo "== $m: guard / artefact =="
    grep -E "GUARD FAILED|ARTEFACT FAILURES|fit file missing" "$L/dr_fits_$m.log" 2>/dev/null | head -8 || echo "  (none)"
  done
  echo; echo "== inclusive plateaus =="
  for s in pp pbpb; do for st in 3 4; do for sub in "" "pt4bin/"; do
    f=~/usatlasdata/dimuon_data/plots/${s}_trigger_efficiency/mc_based/step${st}_dr_fit/${sub}plateau_guard_report.txt
    [ -f "$f" ] && printf "%-6s step%s %-8s %s | %s | %s\n" "$s" "$st" "${sub:-8bin}" \
      "$(grep -m1 '^# inclusive' "$f" | sed 's/# inclusive plateau = //')" \
      "$(grep -m1 '^FLAGGED' "$f" | sed 's/.*: //') flagged" \
      "$(grep -m1 '^FAILING' "$f" | sed 's/.*: //') failing"
  done; done; done
  echo; echo "== PNG counts =="
  echo "  nominal 8-bin: $(find ~/usatlasdata/dimuon_data/plots/{pp,pbpb}_trigger_efficiency/mc_based -name '*.png' -not -path '*pt4bin*' 2>/dev/null | wc -l)"
  echo "  4-bin variant: $(find ~/usatlasdata/dimuon_data/plots/{pp,pbpb}_trigger_efficiency/mc_based -name '*.png' -path '*pt4bin*' 2>/dev/null | wc -l)"
} > "$OUT" 2>&1
