#!/usr/bin/env bash
# One-off (2026-10-13): data NTuple rerun for ONE dataset (pp24 | pbpb23 | pbpb25 | pbpb26), tracking doc
# docs/tracking/ntuple_rerun_pp24_pbpb2526_20261013.md: nominal (crossx) + trigger-efficiency variants,
# Condor submit -> wait (held / >2 h running = FAIL, no silent hang) -> validate -> hadd.
# The previous combined files are copied to <data dir>/pre_rerun_20261013/ first.
# VERIFY_ONLY=<cl1>,<cl2>: skip submission, validate + hadd those already-finished clusters.
# Last line of the log is always "ALL DONE <ds>" or "FATAL <ds>: ...".
source /usatlas/u/yuhanguo/setup.sh >/dev/null 2>&1
set -uo pipefail
DS=${1:?usage: $0 pp24|pbpb23|pbpb25|pbpb26}
A=/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis
if [[ -n ${VERIFY_ONLY:-} && -z ${T0_OVERRIDE:-} ]]; then echo "VERIFY_ONLY needs T0_OVERRIDE (epoch <= submission time)"; exit 1; fi
T0=${T0_OVERRIDE:-$(date +%s)}
BCHK=$A/pipelines/logs/build_check_runlb_${DS}_$$; mkdir -p $BCHK
log(){ echo "[$(date '+%F %T')] $*"; }
fail(){ log "FATAL $DS: $*"; exit 1; }
case $DS in
  pp24)   D=~/usatlasdata/dimuon_data/pp_2024;   NQ=13; PFX=pp_2024
          SUBS=(run_pp_24.sub run_pp_24_nominal.sub)
          TAGS=(single_mu4_mindR_0_02_res_cut_v2 2mu4_mindR_0_02) ;;
  pbpb23) D=~/usatlasdata/dimuon_data/pbpb_2023; NQ=6;  PFX=pbpb_2023   # part6 = re-skimmed lost events of 463315/462759
          SUBS=(run_pbpb_23.sub run_pbpb_23_nominal.sub)
          TAGS=(single_mu4_mindR_0_02_res_cut_v2 single_mu4_mindR_0_02) ;;
  pbpb25) D=~/usatlasdata/dimuon_data/pbpb_2025; NQ=7;  PFX=pbpb_2025
          SUBS=(run_pbpb_25.sub run_pbpb_25_nominal.sub)
          TAGS=(single_mu4_mindR_0_02_res_cut_v2 single_mu4_mindR_0_02) ;;
  pbpb26) D=~/usatlasdata/dimuon_data/pbpb_2026; NQ=7;  PFX=pbpb_2026
          SUBS=(run_pbpb_26.sub run_pbpb_26_nominal.sub)
          TAGS=(single_mu4_mindR_0_02_res_cut_v2 single_mu4_mindR_0_02) ;;
  *) fail "unknown dataset" ;;
esac
for s in "${SUBS[@]}"; do
  q=$(sed -n 's/^queue \([0-9]\+\)$/\1/p' $A/NTupleProcessingCode/$s)
  [[ "$q" == "$NQ" ]] || fail "$s queue=$q != part count $NQ"
done

if [[ -z ${VERIFY_ONLY:-} ]]; then
  # Pb+Pb: the event-selection cuts file must carry the cut-2 window keys (PassEventSel throws otherwise)
  if [[ $DS == pbpb* ]]; then
    f=$D/event_sel_cuts_${PFX}.root
    root -l -b -q -e "auto f=TFile::Open(\"$f\"); if(!f||f->IsZombie()||!f->Get(\"ZDC_time_A_lo_ns\")) gSystem->Exit(3); gSystem->Exit(0);" >/dev/null 2>&1 \
      || fail "$f lacks the cut-2 window keys (stale cuts file)"
    # CUTS_AFTER (epoch): the cuts file must be newer (e.g. the 2026 D8 re-derivation)
    [[ -z ${CUTS_AFTER:-} || $(stat -c %Y $f) -gt $CUTS_AFTER ]] || fail "$f older than CUTS_AFTER"
  fi
  mkdir -p $D/pre_rerun_20261013
  for t in "${TAGS[@]}"; do
    for f in muon_pairs_${PFX}_${t}.root hists_cut_acceptance_${PFX}_${t}.root; do
      [[ -f $D/$f && ! -f $D/pre_rerun_20261013/$f ]] && cp -p $D/$f $D/pre_rerun_20261013/
    done
  done
  log "backup of combined files in $D/pre_rerun_20261013"
  cd $A/NTupleProcessingCode; mkdir -p logs
  CIDS=()
  for s in "${SUBS[@]}"; do
    out="$(condor_submit -append 'periodic_hold = (JobStatus == 2) && ((time() - EnteredCurrentStatus) > 7200)' $s)" || fail "condor_submit $s"
    c="$(echo "$out" | sed -n 's/.*cluster \([0-9]\+\).*/\1/p' | tail -n1)"
    [[ -n $c ]] || fail "no cluster id for $s"
    CIDS+=($c); log "submitted $s cluster $c"
  done
else
  IFS=, read -ra CIDS <<< "$VERIFY_ONLY"
  log "VERIFY_ONLY clusters ${CIDS[*]} (T0=$T0)"
fi
cd $A/NTupleProcessingCode

# wait: a failing condor_q (schedd timeout) is retried, never read as "finished"
for c in "${CIDS[@]}"; do
  nerr=0
  while :; do
    st=$(condor_q $c -autoformat JobStatus 2>/dev/null); rc=$?
    if (( rc != 0 )); then (( ++nerr > 10 )) && fail "condor_q failing for cluster $c"; sleep 60; continue; fi
    nerr=0
    [[ -z $st ]] && break
    held=$(grep -c '^5$' <<< "$st")
    (( held > 0 )) && fail "cluster $c has $held HELD jobs (incl. >2 h running watchdog)"
    sleep 60
  done
  # every proc must have COMPLETED (JobStatus 4) with ExitCode 0
  for try in 1 2 3 4 5; do   # history records can lag the queue
    hist=$(condor_history $c -af ProcId JobStatus ExitCode 2>/dev/null)
    nok=$(awk '$2==4 && $3==0' <<< "$hist" | sort -u -k1,1n | wc -l)
    (( nok == NQ )) && break; sleep 30
  done
  (( nok == NQ )) || fail "cluster $c: only $nok/$NQ procs completed with ExitCode 0: $(tr '\n' ';' <<< "$hist")"
  log "cluster $c finished: $nok/$NQ procs ExitCode 0"
done

# validate: .err only the known TStreamerInfo warnings; every part fresh and readable by content
for s in "${SUBS[@]}"; do
  b=${s%.sub}
  for p in $(seq 0 $((NQ-1))); do
    e=logs/${b}_batch_${p}.err
    bad=$(grep -v 'Warning in <TStreamerInfo::Build>' $e 2>/dev/null | grep -iE 'error|exception|terminate|fatal|segmentation' | head -3)
    [[ -z $bad ]] || fail "$e: $bad"
  done
done
for t in "${TAGS[@]}"; do
  mp=(); h=()
  for p in $(seq 1 $NQ); do
    f=$D/muon_pairs_${PFX}_part${p}_${t}.root; g=$D/hists_cut_acceptance_${PFX}_part${p}_${t}.root
    for x in $f $g; do
      [[ -s $x && $(stat -c %Y $x) -gt $T0 ]] || fail "batch output missing/stale: $x"
    done
    r=$(root -l -b -q -e "auto f=TFile::Open(\"$f\"); auto g=TFile::Open(\"$g\"); if(!f||f->IsZombie()||!g||g->IsZombie()){printf(\"CHK bad\n\");return;} auto t=(TTree*)f->Get(\"muon_pair_tree_sign2\"); auto h=(TH1*)g->Get(\"h_cutAcceptance_sign2\"); printf(\"CHK %lld %g\n\", t? t->GetEntries():-1LL, h? h->GetEntries():-1.);" 2>/dev/null | sed -n 's/^CHK //p')
    read -r nt nh <<< "$r"
    [[ $nt =~ ^[0-9]+$ && $nt -gt 0 ]] && awk "BEGIN{exit !($nh>0)}" || fail "part content check failed for part $p $t: '$r'"
    mp+=($f); h+=($g)
  done
  hadd -f $D/muon_pairs_${PFX}_${t}.root "${mp[@]}" >/dev/null || fail "hadd mp $t"
  hadd -f $D/hists_cut_acceptance_${PFX}_${t}.root "${h[@]}" >/dev/null || fail "hadd h $t"
  n=$(root -l -b -q -e "auto f=TFile::Open(\"$D/muon_pairs_${PFX}_${t}.root\"); auto t=(TTree*)f->Get(\"muon_pair_tree_sign2\"); printf(\"NOS=%lld\n\", t? t->GetEntries() : -1LL);" 2>/dev/null | sed -n 's/NOS=//p')
  [[ -n $n && $n -gt 0 ]] || fail "hadd output $t has no OS pairs ($n)"
  # run/LB exclusion verified on the OUTPUT (§3a): lb filled, no pair in an excluded run/LB
  isp=true; yr=${PFX##*_}; [[ $DS == pp24 ]] && isp=false
  ex=$(cd $A/pipelines && root -l -b -q -e "gSystem->SetBuildDir(\"$BCHK\",true); int e=0; gROOT->ProcessLine(\".L check_runlb_exclusion.C+\",&e); if(e) gSystem->Exit(2); gROOT->ProcessLine(\"check_runlb_exclusion(\\\"$D/muon_pairs_${PFX}_${t}.root\\\", $isp, $yr)\");" 2>&1 < /dev/null); rc=$?
  ex=$(sed -n 's/^EXCL //p' <<< "$ex")
  (( rc == 0 )) || fail "run/LB exclusion check failed for $t (rc=$rc): nlb0 nbad maxlb nruns = '$ex'"
  log "hadd $t OK, OS pairs = $n; exclusion check (nlb0 nbad maxlb nruns) = $ex"
done
log "ALL DONE $DS"
