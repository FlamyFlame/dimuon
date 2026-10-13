#!/bin/bash
# Pb+Pb data-QA refresh for one year (from ntuple_rerun_pp24_pbpb2526_20261013.md step 6; used by run_pbpb23_rerun_part6_20261013.sh): per-run N/L (nominal NTuple) +
# factored N_HLT/L, N_evsel/N_HLT, N_sig/N_evsel (raw NTUP, current cuts file, /pnfs over xrootd).
# usage: qa_pbpb_year.sh <yr> <nparts> <good_lb_only_lumi true|false> <work dir>
source /usatlas/u/yuhanguo/setup.sh >/dev/null 2>&1
export X509_USER_PROXY=/tmp/x509up_u101379
Y=$1; N=$2; GL=$3; S=${4:?work dir}
A=/gpfs/mnt/atlasgpfs01/usatlas/workarea/yuhanguo/dimuon_codes/Analysis; D=$HOME/usatlasdata/dimuon_data; P=$D/plots/data_qa
L=/gpfs/mnt/atlasgpfs01/usatlas/workarea/yuhanguo/dimuon_codes/IntNotes/data/luminosity
B=$S/bqa_pbpb$Y; E=$S/evsel$Y; mkdir -p $B $E; cd $A
fail(){ echo "FATAL pbpb$Y QA: $*"; exit 1; }
dl(){ local t=$(date +%s); "$@" < /dev/null & local p=$!; while kill -0 $p 2>/dev/null; do (( $(date +%s)-t > 10800 )) && { kill $p; fail "deadline: $1"; }; sleep 15; done; wait $p || fail "rc: $*"; }
INC="gSystem->SetBuildDir(\"$B\",true); gSystem->AddIncludePath(\"-I$A\"); gInterpreter->AddIncludePath(\"$A\");"
dl root -l -b -q -e "$INC" -e "int e=0; gROOT->ProcessLine(\".L plotting_codes/data_qa/plot_signal_yield_per_lumi_per_run.cxx+\",&e); if(e) gSystem->Exit(2); gROOT->ProcessLine(\"plot_signal_yield_per_lumi_per_run(\\\"pbpb$Y\\\")\",&e); gSystem->Exit(e?3:0);"
echo "OK N/L"
root -l -b -q -e "$INC" -e "int e=0; gROOT->ProcessLine(\".L plotting_codes/data_qa/evsel_single_run.cxx+\",&e); gSystem->Exit(e?2:0);" > $E/compile.log 2>&1 || fail "compile evsel_single_run"
pids=()
for p in $(seq 1 $N); do
  f=$D/pbpb_20$Y/data_pbpb${Y}_part$p.root; t=$(readlink $f); [[ $t == /pnfs/* ]] && f="root://dcgftp.usatlas.bnl.gov:1096/$t"
  ( dl root -l -b -q -e "$INC" -e ".L plotting_codes/data_qa/evsel_single_run.cxx+" -e "evsel_single_run_fill($Y,\"$f\",\"$E/part$p.root\",{})" > $E/part$p.log 2>&1 ) &
  pids+=($!)
done
for q in "${pids[@]}"; do wait $q || fail "evsel fill part"; done
for p in $(seq 1 $N); do grep -q "^wrote" $E/part$p.log || fail "part$p: $(tail -2 $E/part$p.log)"; done
hadd -f $E/merged.root $(for p in $(seq 1 $N); do echo $E/part$p.root; done) >/dev/null || fail hadd
echo "OK evsel fill"
dl root -l -b -q -e "$INC" -e "int e=0; gROOT->ProcessLine(\".L plotting_codes/data_qa/run_norm_diagnostics.cxx+\",&e); if(e) gSystem->Exit(2); gROOT->ProcessLine(\"run_norm_diagnostics($Y,\\\"$E/merged.root\\\",\\\"$D/pbpb_20$Y/muon_pairs_pbpb_20${Y}_single_mu4_mindR_0_02.root\\\",\\\"$L/pbpb_20$Y/lumitable_pbpb_${Y}_HLT_mu4.csv\\\",{},\\\"$P\\\",10,$GL)\",&e); gSystem->Exit(e?3:0);"
echo "ALL DONE pbpb$Y QA"
