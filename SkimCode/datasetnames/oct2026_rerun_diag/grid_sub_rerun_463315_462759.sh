#!/bin/bash
# DIAGNOSTIC re-skim (2026-10-12) of the two Pb+Pb 2023 runs whose skimmed HLT_mu4 events per
# lumi are low (463315: 0.43x, 462759: 0.74x; tracking doc data_per_run_signal_yield_qa.md R2).
# Same input datasets (r16069_p6447), GRL and run mode as the May-2026 parts 2/3 (tasks
# 50267236/50267250); source = current HEAD (ZDC-skip policy fe04e1f, identical on complete-ZDC events).
# Output name deliberately has NO "partN": it must never be auto-merged into data_pbpb23_partN
# (that would double-count these runs). Compare per-LB counts with parts 2/3 by hand.
# Submit from a shell that sourced ../setup_25.sh (AthAnalysis 25.2.89) + lsetup panda.
pathena --trf "TRIGRATES_RUNMODE=hi2023 athena.py TrigRates_CA.py --filesInput=%IN --evtMax=-1" \
        --inDsTxt InDstxt_PbPb2023data_rerun_463315_462759.txt \
        --outDS user.yuhang.TrigRates.dimuon.PbPb2023data.Oct2026.v1.rerun463315_462759. \
        --mergeOutput \
        --nGBPerJob MAX \
        --extFile data23_hi.periodAllYear_DetStatus-v120-pro33-03_MERGED_PHYS_HeavyIon_All_Good_IgnoreBSPOT_INVALID.xml \
        --extOutFile myfile.root \
        --excludeFile=clean.sh \
        --excludeFile=vim_backup
