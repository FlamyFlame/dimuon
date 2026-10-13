#!/usr/bin/env bash
# One-off (2026-10-12): after the 2026 bad-run / bad-LB exclusion (PbPbEventSelConfig.h,
# tracking doc pbpb2026_analysis_support.md Step 17) rerun for PbPb 2026 ONLY:
#   Stage 0  event-selection derivation + plots (event_level, cuts, cuts_alt) and the
#            cross-year FCal comparison plots;
#   Stage 1  NTuple processing: run_pbpb_26.sub (trig-eff, _res_cut_v2) and
#            run_pbpb_26_nominal.sub (crossx, _mindR_0_02); wait; validate; hadd.
# Deliberately STOPS before RDF / trig-eff fits / crossx (user: do not run crossx yet).
source /usatlas/u/yuhanguo/setup.sh >/dev/null 2>&1
set -uo pipefail
A=/usatlas/u/yuhanguo/workarea/dimuon_codes/Analysis
D=/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pbpb_2026
log(){ echo "[$(date '+%F %T')] $*"; }
fail(){ log "FATAL: $*"; exit 1; }
NQ=7

cd $A/plotting_codes/event_selection
log "Stage 0: event_level(26)"
root -l -b -q "plot_pbpb_event_sel_event_level.cxx+(26)" || fail "event_level rc=$?"
log "Stage 0: cuts(26)"
root -l -b -q "plot_pbpb_event_sel_cuts.cxx+(26)" || fail "cuts rc=$?"
log "Stage 0: cuts_alt(26)"
root -l -b <<ROOTEOF || fail "cuts_alt rc=$?"
.L plot_pbpb_event_sel_cuts.cxx+
plot_pbpb_event_sel_cuts_alt(26);
gSystem->Exit(0);
ROOTEOF
for f in event_sel_cuts_pbpb_2026.root event_sel_cuts_pbpb_2026_alt.root; do
  [[ $D/$f -nt $0 ]] || fail "$f not refreshed"
done
for fn in plot_pbpb_fcal_comparison plot_pbpb_fcal_comparison_alt plot_pbpb_fcal_comparison_262524 plot_pbpb_fcal_comparison_262524_alt; do
  log "Stage 0: $fn"
  root -l -b <<ROOTEOF || fail "$fn rc=$?"
.L plot_pbpb_fcal_comparison.cxx+
$fn();
gSystem->Exit(0);
ROOTEOF
done
log "Stage 0 DONE"

cd $A/NTupleProcessingCode; mkdir -p logs
declare -A CID
for sub in run_pbpb_26.sub run_pbpb_26_nominal.sub; do
  out="$(condor_submit $sub)" || fail "condor_submit $sub"
  CID[$sub]="$(echo "$out" | sed -n 's/.*cluster \([0-9]\+\).*/\1/p' | tail -n1)"
  log "submitted $sub cluster ${CID[$sub]}"
done
for sub in "${!CID[@]}"; do
  while [[ -n "$(condor_q ${CID[$sub]} -autoformat ClusterId 2>/dev/null)" ]]; do
    held=$(condor_q ${CID[$sub]} -autoformat JobStatus 2>/dev/null | grep -c '^5$')
    [[ $held -gt 0 ]] && fail "cluster ${CID[$sub]} has $held HELD jobs"
    sleep 120
  done
  log "cluster ${CID[$sub]} ($sub) finished"
done

for tag in single_mu4_mindR_0_02_res_cut_v2 single_mu4_mindR_0_02; do
  mp=(); h=()
  for p in $(seq 1 $NQ); do
    f=$D/muon_pairs_pbpb_2026_part${p}_${tag}.root; g=$D/hists_cut_acceptance_pbpb_2026_part${p}_${tag}.root
    [[ -s $f && $f -nt $0 && -s $g && $g -nt $0 ]] || fail "batch output missing/stale: $f / $g"
    mp+=($f); h+=($g)
  done
  hadd -f $D/muon_pairs_pbpb_2026_${tag}.root "${mp[@]}" >/dev/null || fail "hadd mp $tag"
  hadd -f $D/hists_cut_acceptance_pbpb_2026_${tag}.root "${h[@]}" >/dev/null || fail "hadd h $tag"
  log "hadd $tag OK"
done
log "ALL DONE"
