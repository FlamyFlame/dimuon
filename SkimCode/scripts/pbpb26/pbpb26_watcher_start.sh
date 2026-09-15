#!/bin/bash
# (Re)start the standalone PbPb2026 watcher exactly once.
#  - kills any previous watcher by SESSION id (the watcher runs under setsid, so its whole
#    session -- main loop + transient subshells -- goes together); never by name pattern.
#  - records the real child PID: `$!` after `setsid nohup ... &` is the short-lived setsid
#    wrapper, not the watcher, so the pid is found via the new session instead.
D=/usatlas/u/yuhanguo/usatlasdata/dimuon_data
PIDF=$D/pbpb26_watcher.pid
old=$(cat "$PIDF" 2>/dev/null)
if [[ -n "$old" ]] && kill -0 "$old" 2>/dev/null; then
  sid=$(ps -o sid= -p "$old" | tr -d ' ')
  [[ -n "$sid" ]] && kill -- -"$sid" 2>/dev/null   # negative pid = whole process group
  for _ in 1 2 3 4 5; do kill -0 "$old" 2>/dev/null || break; sleep 1; done
fi
setsid nohup bash "$D/pbpb26_watcher.sh" > /dev/null 2>&1 &
new=""
for _ in 1 2 3 4 5 6 7 8 9 10; do
  new=$(pgrep -f "bash $D/pbpb26_watcher.sh" | while read -r p; do [[ "$(ps -o ppid= -p $p | tr -d ' ')" == "1" ]] && echo $p; done | head -1)
  [[ -n "$new" ]] && break; sleep 1
done
echo "$new" > "$PIDF"
echo "watcher started pid=$new sid=$(ps -o sid= -p $new | tr -d ' ')"
