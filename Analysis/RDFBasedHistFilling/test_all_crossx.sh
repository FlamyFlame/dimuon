#!/bin/bash
set -e

source ~/setup.sh > /dev/null 2>&1

echo "========== Testing pp24 =========="
root -b -q << 'ROOTEOF'
{
    gSystem->cd("/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/RDFBasedHistFilling");
    .L RDFBasedHistFillingPP.cxx
    RDFBasedHistFillingPP filler_pp24(24, 2);  // run_year=24, trigger_mode=2
    filler_pp24.Run();
}
ROOTEOF

echo "========== Testing pbpb23 =========="
root -b -q << 'ROOTEOF'
{
    gSystem->cd("/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/RDFBasedHistFilling");
    .L RDFBasedHistFillingPbPb.cxx
    RDFBasedHistFillingPbPb filler_pbpb23(23, 2);  // run_year=23, trigger_mode=2
    filler_pbpb23.Run();
}
ROOTEOF

echo "========== Testing pbpb24 =========="
root -b -q << 'ROOTEOF'
{
    gSystem->cd("/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/RDFBasedHistFilling");
    .L RDFBasedHistFillingPbPb.cxx
    RDFBasedHistFillingPbPb filler_pbpb24(24, 2);  // run_year=24, trigger_mode=2
    filler_pbpb24.Run();
}
ROOTEOF

echo "========== All tests completed =========="
