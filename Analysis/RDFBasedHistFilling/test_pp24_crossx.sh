#!/bin/bash
# Test pp24 crossx histogram filling

source ~/setup.sh > /dev/null 2>&1

cd /usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/RDFBasedHistFilling

echo "[TEST] Running pp24 with trigger_mode=2..."
timeout 300 root -b -q << 'EOF' 2>&1
{
    gROOT->SetBatch(kTRUE);
    gROOT->LoadMacro("RDFBasedHistFillingPP.cxx++O"); // Compile and load
    // Wait for compilation to finish
    gSystem->Sleep(1000);
    
    RDFBasedHistFillingPP pp24(24, false);
    pp24.trigger_mode = 2;
    pp24.FillHistograms();
    
    std::cout << "\n[SUCCESS] pp24 completed!" << std::endl;
}
.q
EOF
