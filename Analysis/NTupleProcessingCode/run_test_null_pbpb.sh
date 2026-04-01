#!/bin/bash
cd /usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/NTupleProcessingCode

source ~/setup.sh

root -b -l << EOF
.L DataAnalysisClasses.h

PbPbAnalysis pbpb_23 (23, 1);
pbpb_23.nevents_max = 200;
pbpb_23.trigger_mode = 1;
pbpb_23.Run();

.q
EOF
