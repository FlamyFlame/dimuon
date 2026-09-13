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
# Also cross-checks the FIVE places that must agree on the year's part count:
#   1. file_batch_max{<yr>} in NTupleProcessingCode/PbPbExtras.c
#   2. `queue N` in every run_pbpb_<yr>*.sub (bar the deliberately-partial one-offs)
#   3. QUEUE_COUNTS[<yr>] in both Pb+Pb pipelines
#   4. ScrambGen::NParts(<yr>)  -- builds the mixed-event template T_mix
#   5. the per-part file lists in the six event-selection / preamp / FCal plotting macros
# and, when the year's directory exists, all of them against the files actually on disk.
#
# The two directions are NOT equally dangerous, and classes 2 and 5 are graded accordingly:
# an UNDER-count silently processes a subset of the data (error), an OVER-count throws
# loudly on the missing part (warning).  Classes 1, 3, 4 and the on-disk comparison are
# exact-match -- for those, either direction is silent.
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
if [[ "${nmax:-0}" == "0" ]]; then
  # A deliberate "unset" (see PbPbExtras.c): the year's skim is still in production, so its
  # part count is not knowable yet.  Cross-checking the other declarations against 0 would
  # emit several mismatches that all say the same thing, so say it once instead.
  echo "  >> Pb+Pb 20${YR} is deliberately UNSET: its part count is not established yet, so the"
  echo "     other declarations (queue N, QUEUE_COUNTS, ScrambGen::NParts, plotting part lists)"
  echo "     are moot until it is. Set them TOGETHER from the files on disk, then re-run this."
  bad=1
else
bad=0
# Every run_pbpb_<yr>*.sub that processes the year's FULL file set must agree: a queue count
# below the part count SILENTLY processes only part of the data.  The deliberately-partial
# one-offs are named explicitly and excluded, so a NEW variant is checked by default rather
# than quietly exempt.  (grid_monitor auto-updates only run_pbpb_<yr>.sub, so the others are
# exactly the ones that drift.)
partial_variants=( "_part4_only" "_test" "_diagnose" )
for f in "$A"/NTupleProcessingCode/run_pbpb_${YR}*.sub; do
  [[ -f "$f" ]] || continue
  base=$(basename "$f"); skip=0
  for pv in "${partial_variants[@]}"; do [[ "$base" == *"${pv}"* ]] && skip=1; done
  if [[ $skip -eq 1 ]]; then
    q=$(grep -oP '^queue \K[0-9]+' "$f" | head -1 || true)
    echo "  note     ${base}: queue $q (deliberately partial variant, not checked)"
    continue
  fi
  q=$(grep -oP '^queue \K[0-9]+' "$f" | head -1 || true)
  # UNDER-count is the dangerous direction: Condor silently processes only part of the data.
  # OVER-count throws loudly on the missing part, so it is stale-but-safe -> warn, don't fail.
  if [[ "$q" -lt "$nmax" ]]; then
    echo "  MISMATCH ${base}: queue $q < $nmax -- SILENTLY drops data"; bad=1
  elif [[ "$q" -gt "$nmax" ]]; then
    echo "  warn     ${base}: queue $q > $nmax (stale; throws loudly on the missing part)"
  fi
done

# ScrambGen::NParts(<yr>) -- builds the mixed-event template T_mix from the single-muon trees.
sg=$(grep -oP "case ${YR}: return \K[0-9]+" "$A/ScrambGen/ScrambGen.h" | head -1 || true)
if [[ -n "$sg" && "$sg" != "$nmax" ]]; then
  echo "  MISMATCH ScrambGen::NParts(${YR}) = $sg != $nmax"; bad=1
elif [[ -z "$sg" ]]; then
  echo "  MISMATCH ScrambGen::NParts has no case for year ${YR}"; bad=1
fi

# Per-part file lists in the event-selection / preamp / FCal plotting macros.  Each hard-codes
# data_pbpb<yr>_part1..N.root; a list shorter than the part count silently derives that year's
# cuts from a subset of its data.
for f in "$A"/plotting_codes/event_selection/plot_pbpb_event_sel_cuts.cxx \
         "$A"/plotting_codes/event_selection/plot_pbpb_event_sel_event_level.cxx \
         "$A"/plotting_codes/event_selection/plot_pbpb_fcal_comparison.cxx \
         "$A"/plotting_codes/event_selection/plot_zdc_preamp_cut_vs_zdcamp.cxx \
         "$A"/plotting_codes/event_selection/plot_zdc_preamp_cut_over_mean.cxx \
         "$A"/plotting_codes/event_selection/plot_zdc_preamp_gauss_fit.cxx; do
  [[ -f "$f" ]] || continue
  # Count DISTINCT path literals, not matching lines: a line count includes comments, so a
  # single comment naming data_pbpb<yr>_partN.root would inflate the total -- and in the
  # dangerous direction it could lift a genuine UNDER-count back to the expected number and
  # silence this very check.  Strip // comments first, then count unique part filenames.
  cnt=$(sed 's://.*::' "$f" | grep -o "data_pbpb${YR}_part[0-9]\+\.root" | sort -u | wc -l)
  if [[ "$cnt" -gt 0 && "$cnt" -lt "$nmax" ]]; then
    echo "  MISMATCH $(basename "$f"): lists $cnt of $nmax part file(s) for 20${YR}"\
         "-- derives from a SUBSET of the data"; bad=1
  elif [[ "$cnt" -gt "$nmax" ]]; then
    echo "  warn     $(basename "$f"): lists $cnt part file(s) for 20${YR}, expected $nmax"
  fi
done
for p in pipeline_pbpb_crossx.sh pipeline_pbpb_trig_eff.sh; do
  qc=$(grep -oP "\[${YR}\]=\K[0-9]+" "$SCRIPT_DIR/$p" | head -1 || true)
  if [[ "$qc" != "$nmax" ]]; then echo "  MISMATCH $p: QUEUE_COUNTS[$YR]=$qc != $nmax"; bad=1; fi
done
if [[ -d "$DIR" ]]; then
  # Match data_pbpb<yr>_part<digits>.root EXACTLY.  A plain glob also matches the
  # data_pbpb<yr>_part1.bak_<date>.root files grid_monitor leaves while re-merging, which
  # would inflate the count -- and an inflated count hides a genuine shortfall.
  mapfile -t realparts < <(ls "$DIR" 2>/dev/null | grep -E "^data_pbpb${YR}_part[0-9]+\\.root$" | sort -V)
  non=${#realparts[@]}
  echo "  part files on disk = $non"
  if [[ $non -gt 0 ]]; then
    # Part numbers must be CONTIGUOUS 1..N: the Condor model submits one job per
    # file_batch in 1..queue, so a hole (e.g. parts 1,4 with 2,3 absent) means those jobs
    # die on a missing input while the count still looks plausible.
    nums=(); for f in "${realparts[@]}"; do n="${f#*part}"; nums+=( "${n%%.root}" ); done
    expect=1; holes=()
    for n in "${nums[@]}"; do
      while [[ $expect -lt $n ]]; do holes+=("$expect"); expect=$((expect+1)); done
      expect=$((n+1))
    done
    if [[ ${#holes[@]} -gt 0 ]]; then
      echo "  MISMATCH part numbering is NOT contiguous: present ${nums[*]}; missing ${holes[*]}."
      echo "           Condor submits one job per batch 1..N, so the missing batches would fail."
      echo "           While a skim is still recovering this can be TRANSIENT (a late part has not"
      echo "           landed yet) -- this check is a pre-submission gate, so the answer is 'wait',"
      echo "           not 'renumber'. Confirm the final set with the skimming session."
      bad=1
    fi
    if [[ -n "${nmax:-}" && "$non" != "$nmax" ]]; then
      echo "  MISMATCH: $non part file(s) on disk but file_batch_max is ${nmax}"; bad=1
    fi
  fi
else
  echo "  (no $DIR yet -- on-disk count not checked)"
fi
[[ $bad -eq 0 ]] && echo "  all declared part counts agree"

fi

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
