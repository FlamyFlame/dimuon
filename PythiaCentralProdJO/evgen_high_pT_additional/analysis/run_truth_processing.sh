#!/bin/bash
# Ntuple processing (PythiaTruthAnalysis, external-truth raw-counts mode; Analysis/docs/pythia_truth.md)
# on the downloaded grid-test NTUPs, one slice per config. Outputs:
#   $OUT/muon_pairs_pythia_truth_<label>_kn<N>_no_data_resonance_cuts.root  (consumed by
#   Analysis/plotting_codes/trig_effcy/mc_based/mc_request_stats_boost.cxx)
TP_BASE=${1:-/usatlas/u/yuhanguo/usatlasdata/evgen_high_pT_additional_test}
OUT=${2:-$TP_BASE/truth_processing}
mkdir -p $OUT
source /usatlas/u/yuhanguo/setup.sh > /dev/null 2>&1
cd /gpfs/mnt/atlasgpfs01/usatlas/workarea/yuhanguo/dimuon_codes/Analysis/NTupleProcessingCode
#    config dir                 kn  label (E_COM in the label: the mode's output name has no ecom tag)
for line in "run_sign_indep_pT70_125  4 5p36_SI_m10_pT45" "run_sign_indep_pT125_300 5 5p36_SI_m10_pT45" \
            "run_same_sign_pT24_40    2 5p36_SS_m10"      "run_same_sign_pT40_70    3 5p36_SS_m10_pT20" \
            "run_same_sign_pT70_125   4 5p36_SS_m10_pT20" "run_same_sign_pT125_300  5 5p36_SS_m10_pT20"; do
  read cfg kn label <<< "$line"
  files=$(find $TP_BASE/$cfg -name '*.root' ! -name '*log*' | sort | paste -sd, -)
  [ -z "$files" ] && { echo "SKIP $cfg: no NTUP"; continue; }
  echo "== $cfg kn=$kn label=$label files: $files"
  root -l -b -q -e ".L PythiaAnalysisClasses.h" -e "PythiaTruthAnalysis py($kn,false,5.36); py.truth_input_override=\"$files\"; py.truth_output_dir_override=\"$OUT\"; py.truth_output_label=\"$label\"; py.truth_raw_counts=true; py.Run();" > $OUT/log_${label}_kn${kn}.txt 2>&1
  echo "   rc=$? log $OUT/log_${label}_kn${kn}.txt"
done
