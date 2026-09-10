#!/usr/bin/env bash
# Pre-flight for a Pb+Pb data-taking year: does every artifact the Pb+Pb pipelines
# reference for that year exist?
#
#   ./preflight_pbpb_year.sh 26
#
# Two classes are reported separately, because they fail for different reasons:
#   MISSING  a repo artifact (condor .sub/.sh, RDF runner) that must exist before the
#            pipeline can even be launched -> exit 1.
#   PENDING  a data artifact produced by an earlier stage.  Expected to be absent
#            before that stage has run; NOT an error.
#
# Also cross-checks the three places that must agree on the year's part count
# (file_batch_max in PbPbExtras.c, `queue N` in every run_pbpb_<yr>*.sub, and
# QUEUE_COUNTS in both pipelines).  A DISAGREEMENT here is the one failure that
# silently processes only part of the data, so it is reported as an error.
set -Eeuo pipefail

YR="${1:-}"
[[ -n "$YR" ]] || { echo "usage: $0 <2-digit year, e.g. 26>"; exit 2; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
A="$(cd "${SCRIPT_DIR}/.." && pwd)"
D="/usatlas/u/yuhanguo/usatlasdata/dimuon_data"
DIR="${D}/pbpb_20${YR}"

ok=0; miss=0; pend=0
chk(){     if [[ -e "$2" ]]; then printf '  OK      %-28s %s\n' "$1" "$2"; ok=$((ok+1));
           else printf '  MISSING %-28s %s\n' "$1" "$2"; miss=$((miss+1)); fi; }
pending(){ if [[ -e "$2" ]]; then printf '  PRESENT %-28s %s\n' "$1" "$2"; ok=$((ok+1));
           else printf '  PENDING %-28s %s\n' "$1" "$2"; pend=$((pend+1)); fi; }

echo "=== Pb+Pb 20${YR} pipeline pre-flight ==="
echo "--- repo artifacts (must exist before launching) ---"
for v in "" _nominal _no_res_cut _output_single_muon_tree _mindR_0_01 _no_res_mu4_mu4noL1; do
  chk "condor submit" "$A/NTupleProcessingCode/run_pbpb_${YR}${v}.sub"
  chk "condor exec"   "$A/NTupleProcessingCode/run_pbpb_${YR}${v}.sh"
done
chk "RDF crossx runner"      "$A/RDFBasedHistFilling/run_crossx_hist_filling_pbpb${YR}.sh"
chk "RDF template fit"       "$A/RDFBasedHistFilling/run_template_fit_pbpb${YR}.sh"
chk "RDF template fit mixed" "$A/RDFBasedHistFilling/run_template_fit_mixed_pbpb${YR}.sh"

echo "--- part-count agreement (a mismatch SILENTLY drops data) ---"
# Read the map ENTRY, not a comment that happens to mention "{23, 3}".  The entry is the
# line that declares run_year_to_file_batch_max_map's initializer.
nmax=$(grep -v '^[[:space:]]*//' "$A/NTupleProcessingCode/PbPbExtras.c" \
       | grep -oP "\{${YR}, \K[0-9]+" | head -1 || true)
echo "  PbPbExtras.c file_batch_max{${YR}} = ${nmax:-<absent>}"
bad=0
# The pipelines launch exactly these two; a wrong queue count in them drops data silently.
# Other run_pbpb_<yr>_*.sub are manual/one-off variants (test jobs, the unmaintained
# mu4_mu4noL1 mode, part-N-only reruns) whose queue count is deliberately different --
# report them as advisory, never as a failure.
for v in "" _nominal; do
  f="$A/NTupleProcessingCode/run_pbpb_${YR}${v}.sub"
  [[ -f "$f" ]] || continue
  q=$(grep -oP '^queue \K[0-9]+' "$f" | head -1 || true)
  if [[ "$q" != "$nmax" ]]; then echo "  MISMATCH $(basename "$f"): queue $q != $nmax"; bad=1; fi
done
for f in "$A"/NTupleProcessingCode/run_pbpb_${YR}_*.sub; do
  [[ -f "$f" ]] || continue
  [[ "$f" == *_nominal.sub ]] && continue
  q=$(grep -oP '^queue \K[0-9]+' "$f" | head -1 || true)
  if [[ "$q" != "$nmax" ]]; then echo "  note     $(basename "$f"): queue $q (pipeline count $nmax)"; fi
done
for p in pipeline_pbpb_crossx.sh pipeline_pbpb_trig_eff.sh; do
  qc=$(grep -oP "\[${YR}\]=\K[0-9]+" "$SCRIPT_DIR/$p" | head -1 || true)
  if [[ "$qc" != "$nmax" ]]; then echo "  MISMATCH $p: QUEUE_COUNTS[$YR]=$qc != $nmax"; bad=1; fi
done
if [[ -d "$DIR" ]]; then
  non=$(ls "$DIR"/data_pbpb${YR}_part*.root 2>/dev/null | wc -l)
  echo "  part files on disk = $non"
  if [[ "$non" -gt 0 && "$non" != "$nmax" ]]; then
    echo "  MISMATCH: $non part files on disk but file_batch_max is $nmax"; bad=1
  fi
else
  echo "  (no $DIR yet -- on-disk count not checked)"
fi
[[ $bad -eq 0 ]] && echo "  all declared part counts agree"

echo "--- data artifacts (PENDING is fine before the producing stage runs) ---"
pending "raw skim dir" "$DIR"
if [[ -n "${nmax:-}" ]]; then
  for ((p=1; p<=nmax; ++p)); do pending "raw skim NTUP part$p" "$DIR/data_pbpb${YR}_part${p}.root"; done
fi
pending "event-sel cuts" "$DIR/event_sel_cuts_pbpb_20${YR}.root"
pending "trig-eff pT fits" "$DIR/trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit.root"
pending "crossx hists" "$DIR/histograms_real_pairs_pbpb_20${YR}_single_mu4_no_trg_plots_nominal.root"

echo
echo "SUMMARY: ok=$ok  missing(repo)=$miss  pending(data)=$pend  count-mismatch=$bad"
[[ $miss -eq 0 && $bad -eq 0 ]]
