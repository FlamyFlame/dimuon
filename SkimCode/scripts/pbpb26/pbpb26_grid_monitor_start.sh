#!/bin/bash
# (Re)start grid_monitor.sh for the PbPb2026 skim and record its PID.
#
# NEVER identify it with `pkill -f "grid_monitor.sh ..."`: that pattern also matches the
# command line of whatever wrapper/monitor is doing the pkill, so the caller kills itself.
# (It did exactly that once here -- the release monitor died with exit 144 and took
# grid_monitor with it.)  Use the recorded PID.
D=/usatlas/u/yuhanguo/usatlasdata/dimuon_data
SK=/usatlas/u/yuhanguo/workarea/dimuon_codes/SkimCode/scripts
PIDF=$D/pbpb26_grid_monitor.pid

old=$(cat "$PIDF" 2>/dev/null)
if [[ -n "$old" ]] && kill -0 "$old" 2>/dev/null; then
  kill "$old" 2>/dev/null
  for _ in $(seq 1 20); do kill -0 "$old" 2>/dev/null || break; sleep 1; done
  echo "stopped previous grid_monitor pid=$old"
fi
cd "$SK" || exit 1
nohup bash grid_monitor.sh -i 20 "$D/sep2026_pbpb26_skim.txt" >> "$D/grid_monitor_pbpb26_nohup.log" 2>&1 &
echo $! > "$PIDF"
echo "grid_monitor started pid=$(cat "$PIDF") tasks=$(grep -vcE '^\s*#|^\s*$' "$D/sep2026_pbpb26_skim.txt")"
