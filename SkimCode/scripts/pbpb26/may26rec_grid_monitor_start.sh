#!/bin/bash
D=/usatlas/u/yuhanguo/usatlasdata/dimuon_data; SK=/usatlas/u/yuhanguo/workarea/dimuon_codes/SkimCode/scripts
PIDF=$D/may26rec_grid_monitor.pid
old=$(cat "$PIDF" 2>/dev/null); [[ -n "$old" ]] && kill -0 "$old" 2>/dev/null && { kill "$old"; sleep 2; }
cd "$SK" && setsid nohup bash grid_monitor.sh -i 20 "$D/sep2026_may26_recovery.txt" >> "$D/grid_monitor_may26rec_nohup.log" 2>&1 &
sleep 2; echo $! > "$PIDF"; echo "may26rec grid_monitor started pid=$(cat $PIDF)"
