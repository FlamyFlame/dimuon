#!/bin/bash
# Comprehensive test of FillHistogramsCrossx for all three datasets

source ~/setup.sh > /dev/null 2>&1
cd /usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis/RDFBasedHistFilling

echo "========================================="
echo "Testing FillHistogramsCrossx (trigger_mode==2)"
echo "========================================="

# Test pp24
echo -e "\n[1/3] Testing pp24..."
rm -f histograms_real_pairs_pp_2024_mu4_mu4noL1*.root
mkdir -p logs
timeout 300 root -b -q << 'EOPP24' > logs/pp24_test.log 2>&1
{
  gROOT->SetBatch(kTRUE);
  gROOT->LoadMacro("RDFBasedHistFillingPP.cxx++O");
  spaceData = new RDFBasedHistFillingPP(24, false);
  pp24->trigger_mode = 2;
  pp24->FillHistograms();
  std::cout << "\n[SUCCESS] pp24 completed!" << std::endl;
}
.q
EOPP24
echo "pp24 done. Checking output..."
ls -lh /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/histograms_real_pairs_pp_2024_mu4_mu4noL1*.root 2>&1 | head -5

# Test pbpb23
echo -e "\n[2/3] Testing pbpb23..."
rm -f histograms_real_pairs_pbpb_2023_mu4_mu4noL1*.root
timeout 300 root -b -q << 'EOPBPB23' > logs/pbpb23_test.log 2>&1
{
  gROOT->SetBatch(kTRUE);
  gROOT->LoadMacro("RDFBasedHistFillingPbPb.cxx++O");
  RDFBasedHistFillingPP* pp24 = new RDFBasedHistFillingPP(24, false);
  pbpb23->trigger_mode = 2;
  pbpb23->FillHistograms();
  std::cout << "\n[SUCCESS] pbpb23 completed!" << std::endl;
}
.q
EOPBPB23
echo "pbpb23 done. Checking output..."
ls -lh /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2023/histograms_real_pairs_pbpb_2023_mu4_mu4noL1*.root 2>&1 | head -5

# Test pbpb24
echo -e "\n[3/3] Testing pbpb24..."
rm -f histograms_real_pairs_pbpb_2024_mu4_mu4noL1*.root
timeout 300 root -b -q << 'EOPBPB24' > logs/pbpb24_test.log 2>&1
{
  gROOT->SetBatch(kTRUE);
  gROOT->LoadMacro("RDFBasedHistFillingPbPb.cxx++O");
  RDFBasedHistFillingPbPb* pbpb24 = new RDFBasedHistFillingPbPb(24, false);
  pbpb24->trigger_mode = 2;
  pbpb24->FillHistograms();
  std::cout << "\n[SUCCESS] pbpb24 completed!" << std::endl;
}
.q
EOPBPB24
echo "pbpb24 done. Checking output..."
ls -lh /usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2024/histograms_real_pairs_pbpb_2024_mu4_mu4noL1*.root 2>&1 | head -5

echo -e "\n========================================="
echo "All tests completed!"
echo "========================================="
