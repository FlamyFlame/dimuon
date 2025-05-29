#!/usr/bin/env bash
# ── map pTHatX_YGeV → events_after_filter ─────────────────────────────────────
declare -A afterFilter=(
  [5_8]=1
  [8_14]=10
  [14_24]=50
  [24_40]=200
  [40_70]=500
  [70_125]=2000
  [125_300]=2000
)

# ── job list ──────────────────────────────────────────────────────────────────
jobs=(
user.yuhang.PyJO.dimuon.May25.nn.pTHat125_300GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.np.pTHat125_300GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.pn.pTHat125_300GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.pp.pTHat125_300GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.nn.pTHat70_125GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.np.pTHat70_125GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.pn.pTHat70_125GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.pp.pTHat70_125GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.nn.pTHat40_70GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.np.pTHat40_70GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.pn.pTHat40_70GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.pp.pTHat40_70GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.nn.pTHat24_40GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.nn.pTHat24_40GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.np.pTHat24_40GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.np.pTHat24_40GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.pn.pTHat24_40GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.pn.pTHat24_40GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.pp.pTHat24_40GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.pp.pTHat24_40GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.nn.pTHat14_24GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.nn.pTHat14_24GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.np.pTHat14_24GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.np.pTHat14_24GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.pn.pTHat14_24GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.pn.pTHat14_24GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.pp.pTHat14_24GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.pp.pTHat14_24GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.nn.pTHat8_14GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.nn.pTHat8_14GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.np.pTHat8_14GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.np.pTHat8_14GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.pn.pTHat8_14GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.pn.pTHat8_14GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.pp.pTHat8_14GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.pp.pTHat8_14GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.nn.pTHat5_8GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.nn.pTHat5_8GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.np.pTHat5_8GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.np.pTHat5_8GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.pn.pTHat5_8GeV.job2.v6
user.yuhang.PyJO.dimuon.May25.pn.pTHat5_8GeV.job1.v6
user.yuhang.PyJO.dimuon.May25.pp.pTHat5_8GeV.job1.v6
)

# ── header (optional) ─────────────────────────────────────────────────────────
printf "#X-Y #events_after_filter #events filter_effcy #seconds time/event[s]" > job_metrics.txt
printf "\n" >> job_metrics.txt

for job in "${jobs[@]}"; do
    # extract the pTHat key, e.g. 24_40
    key=$(echo "$job" | grep -oP 'pTHat\K[0-9_]+(?=GeV)')
    after=${afterFilter[$key]:-0}

    # numbers from logs
    events=$(grep -h "events read so far" "${job}_eventLoopHeartBeat.txt"/*txt \
             | tail -1 | grep -o '[0-9]\+' | head -1)
    secs=$(grep -h "Job completed in" "${job}_output.txt"/*txt \
           | tail -1 | awk '{print $(NF-1)}')

    [[ -z $events || -z $secs || $after -eq 0 ]] && \
        { echo "WARN: $job missing data" >&2; continue; }

    tpe=$(awk -v s="$secs" -v e="$events" 'BEGIN{printf("%.3f", s/e)}')
    eff=$(awk -v a="$after" -v e="$events" 'BEGIN{printf("%.3e", a/e)}')

    # replace '_' with '-' for the printed pTHat range
    range=${key//_/-}

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$range" "$after" "$events" "$eff" "$secs" "$tpe" >> job_metrics.txt
done
