#!/usr/bin/env bash
set -Eeuo pipefail

# =============================================================================================
# DeltaR-CORRECTION FIT PIPELINE  (mc_trigger_efficiency.md §3.3 Step 3 / §3.4 Step 4,
# round-7 Autonomy-Contract item 6)
#
# One integrated flow, in the order the physics requires:
#
#   Stage 1  MEASURE + PLATEAU   plot_mc_trig_eff(sample, wp)
#              re-makes the Step-1..4 plots AND writes the per-(pair pT, pair eta) large-dR
#              plateau to  <mc_dir>/dr_correction_plateaus_<label><wp>.root
#   Stage 2  GUARD + FIT         fit_dr_corrections(sample, wp, step, method)
#              reads that ROOT file (never a .txt / .md), enforces |plateau-1| <= 0.1 on a FULL
#              production, divides each cell's curve by ITS OWN plateau, fits
#   Stage 3  PLOT + READ-BACK    plot_dr_correction_fits(sample, wp, step, method)
#              re-opens the fit file in a SEPARATE process, draws measured (black) vs fitted
#              (red) -- 1 PNG per pair-pT bin, 1 subplot per pair-eta bin -- and verifies every
#              persisted function outside its fit range
#
# WHY ARTEFACT VALIDATION AND NOT EXIT CODES: a ROOT macro that throws still exits 0 (the
# exception aborts the interpreter after ROOT has already decided the batch job "ran"). Every
# stage below is therefore checked by looking at what it actually produced.
#
# Usage:
#   ./run_dr_correction_fits.sh
#
# Optional env vars:
#   SAMPLES="pp_full overlay"        samples to process (add `noovl` for the r17663 diagnostic)
#   WPS="tight medium"               working points
#   STEPS="3 4"                      3 = cross / 2mu4 term, 4 = single leg
#   METHODS="..."                    fit methods; one output subdirectory each
#   SKIP_MEASURE=1                   reuse the existing plateau ROOT files (skip Stage 1)
#   SKIP_FIT=1                       reuse the existing fit ROOT files (skip Stage 2) -- for
#                                    re-making only the plots / the read-back audit
#   STRICT_GUARD=1                   a FULL-sample plateau violation stops that sample/step dead
#                                    (no override re-run). Default 0: the violation is reported,
#                                    the fits are re-run with an explicit override so a human can
#                                    LOOK at the plots, and the script still exits non-zero.
# =============================================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
MC_PLOT_DIR="${ANALYSIS_DIR}/plotting_codes/trig_effcy/mc_based"

SAMPLES="${SAMPLES:-pp_full overlay}"
WPS="${WPS:-tight medium}"
STEPS="${STEPS:-3 4}"
METHODS="${METHODS:-powerlaw_fixedRp powerlaw_floatRp expo polyu_fixedRp interp}"
SKIP_MEASURE="${SKIP_MEASURE:-0}"
SKIP_FIT="${SKIP_FIT:-0}"
STRICT_GUARD="${STRICT_GUARD:-0}"

now() { date '+%F %T'; }
log() { echo "[$(now)] $*"; }
fail() { echo "[$(now)] ERROR: $*" >&2; exit 1; }

on_error() {
  local exit_code=$?
  echo "[$(now)] ERROR: command failed (exit=${exit_code}) at line ${BASH_LINENO[0]}: ${BASH_COMMAND}" >&2
}
trap on_error ERR

source_env_once() {
  local err_trap_saved
  err_trap_saved="$(trap -p ERR || true)"
  trap - ERR
  set +e; set +u
  source ~/setup.sh
  local rc=$?
  set -e; set -u
  if [[ -n "${err_trap_saved}" ]]; then eval "${err_trap_saved}"; else trap on_error ERR; fi
  [[ $rc -eq 0 ]] || fail "Failed to source ~/setup.sh"
}

# ---- artefact checks -------------------------------------------------------------------------

# A ROOT file that opens, is non-empty, and contains every named object.
root_file_has_objects() {
  local f="$1"; shift
  [[ -s "$f" ]] || return 1
  local names=("$@") checks=""
  for n in "${names[@]}"; do
    checks+="if (!fin->Get(\"${n}\")) { fin->Close(); gSystem->Exit(4); }"$'\n'
  done
  root -l -b -q <<EOF >/dev/null 2>&1
TFile *fin = TFile::Open("$f", "READ");
if (!fin || fin->IsZombie()) { gSystem->Exit(2); }
${checks}
fin->Close();
gSystem->Exit(0);
EOF
}

wp_flag()   { [[ "$1" == "tight" ]] && echo "true" || echo "false"; }
wp_suffix() { [[ "$1" == "tight" ]] && echo ""     || echo "_medium_wp"; }
wp_dir()    { [[ "$1" == "tight" ]] && echo ""     || echo "medium/"; }

# Sample identity mirrored from dr_correction_sample_cfg.h. Kept minimal (only what the SHELL
# needs to locate artefacts); the macros themselves always read the header, so there is exactly
# one place where a path can be wrong, and it is not this file.
sample_mc_dir() {
  case "$1" in
    pp)      echo "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_test_sample/" ;;
    pp_full) echo "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/" ;;
    overlay) echo "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/" ;;
    noovl)   echo "/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_no_overlay_test_sample/" ;;
    *) fail "unknown sample '$1'" ;;
  esac
}
sample_label() {
  case "$1" in
    pp)      echo "pp24" ;;
    pp_full) echo "pp24_full" ;;
    overlay) echo "hijing_overlay_pbpb23" ;;
    noovl)   echo "r17663_no_overlay" ;;
    *) fail "unknown sample '$1'" ;;
  esac
}
sample_plot_base() {
  case "$1" in
    pp|pp_full) echo "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/mc_based/" ;;
    overlay)    echo "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pbpb_trigger_efficiency/mc_based/" ;;
    noovl)      echo "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/r17663_no_overlay_trigger_efficiency/mc_based/" ;;
    *) fail "unknown sample '$1'" ;;
  esac
}

# ---- go ---------------------------------------------------------------------------------------

log "sourcing ROOT environment"
source_env_once
command -v root >/dev/null 2>&1 || fail "root not on PATH after sourcing ~/setup.sh"

cd "${MC_PLOT_DIR}" || fail "cannot cd to ${MC_PLOT_DIR}"

# ---- Stage 0: compile ONCE ---------------------------------------------------------------------
# `root -q 'X.cxx+'` with no arguments compiles AND RUNS X() with its defaults; `-e '.L X.cxx+'`
# only compiles. Doing it here means the per-sample loops below just dlopen the .so.
log "Stage 0: compiling the three macros with ACLiC (once)"
for m in plot_mc_trig_eff fit_dr_corrections plot_dr_correction_fits; do
  root -l -b -q -e ".L ${m}.cxx+" >/dev/null 2>&1 || true
  [[ -s "${m}_cxx.so" ]] || fail "ACLiC did not produce ${m}_cxx.so -- compile it by hand to see the errors"
  log "  compiled ${m}_cxx.so"
done

declare -a GUARD_FAILURES=()
declare -a ARTEFACT_FAILURES=()
declare -a CHI2_SUMMARY=()

for sample in ${SAMPLES}; do
  MCDIR="$(sample_mc_dir "${sample}")"
  LABEL="$(sample_label "${sample}")"
  PLOTBASE="$(sample_plot_base "${sample}")"

  for wp in ${WPS}; do
    WPF="$(wp_flag "${wp}")"
    WPS_SUF="$(wp_suffix "${wp}")"
    WPD="$(wp_dir "${wp}")"
    PLATEAU_FILE="${MCDIR}dr_correction_plateaus_${LABEL}${WPS_SUF}.root"

    # ---- Stage 1: measure + write the plateau ROOT file ---------------------------------------
    if [[ "${SKIP_MEASURE}" == "1" ]]; then
      log "Stage 1 [${sample}/${wp}]: SKIPPED (SKIP_MEASURE=1), reusing ${PLATEAU_FILE}"
    else
      log "Stage 1 [${sample}/${wp}]: plot_mc_trig_eff -> Step 1-4 plots + plateau ROOT file"
      root -l -b -q "plot_mc_trig_eff.cxx+(\"${sample}\", ${WPF})" > \
        "/tmp/drfit_measure_${sample}_${wp}.log" 2>&1 || true
    fi
    if ! root_file_has_objects "${PLATEAU_FILE}" h_step3_plateau h_step3_plateau_inclusive \
                                                 h_step4_plateau h_step4_plateau_inclusive; then
      ARTEFACT_FAILURES+=("plateau file missing/incomplete: ${PLATEAU_FILE}")
      log "  !! plateau file missing or incomplete -- skipping ${sample}/${wp}"
      continue
    fi
    log "  plateau file OK: ${PLATEAU_FILE}"

    for step in ${STEPS}; do
      GUARD_REPORT="${PLOTBASE}step${step}_dr_fit/${WPD}plateau_guard_report.txt"

      for method in ${METHODS}; do
        FIT_FILE="${MCDIR}dr_correction_fits_${LABEL}${WPS_SUF}_step${step}_${method}.root"
        MDIR="${PLOTBASE}step${step}_dr_fit/${method}/${WPD}"

        # ---- Stage 2: guard + fit --------------------------------------------------------------
        if [[ "${SKIP_FIT}" == "1" ]]; then
          log "Stage 2 [${sample}/${wp}/step${step}/${method}]: SKIPPED (SKIP_FIT=1)"
        else
          # deleted first so the artefact check below cannot pass on a stale file
          rm -f "${FIT_FILE}"
          log "Stage 2 [${sample}/${wp}/step${step}/${method}]: guard + fit"
          root -l -b -q "fit_dr_corrections.cxx+(\"${sample}\", ${WPF}, ${step}, \"${method}\", false)" \
            > "/tmp/drfit_fit_${sample}_${wp}_${step}_${method}.log" 2>&1 || true
        fi

        # The guard writes its verdict BEFORE it throws, so the verdict is readable either way.
        GUARD_FAILED=0
        if [[ -f "${GUARD_REPORT}" ]] && grep -q "^verdict: FAIL" "${GUARD_REPORT}"; then
          GUARD_FAILED=1
        fi

        if [[ "${GUARD_FAILED}" == "1" ]]; then
          if [[ ! " ${GUARD_FAILURES[*]-} " == *" ${sample}/${wp}/step${step} "* ]]; then
            GUARD_FAILURES+=("${sample}/${wp}/step${step}")
          fi
          echo "==============================================================================="
          echo " PLATEAU GUARD FAILED: ${sample} (FULL production) / ${wp} / step ${step}"
          sed -n '/^cells with/,$p' "${GUARD_REPORT}"
          echo "==============================================================================="
          if [[ "${STRICT_GUARD}" == "1" ]]; then
            log "  STRICT_GUARD=1 -> no fits produced for ${sample}/${wp}/step${step}"
            continue
          fi
          if [[ "${SKIP_FIT}" != "1" ]]; then
            log "  re-running WITH the explicit override so the plots exist for human inspection"
            root -l -b -q "fit_dr_corrections.cxx+(\"${sample}\", ${WPF}, ${step}, \"${method}\", true)" \
              > "/tmp/drfit_fit_${sample}_${wp}_${step}_${method}.log" 2>&1 || true
          fi
        fi

        if ! root_file_has_objects "${FIT_FILE}" "h_step${step}_chi2ndf" \
                                                 "h_step${step}_plateau"; then
          ARTEFACT_FAILURES+=("fit file missing/incomplete: ${FIT_FILE}")
          log "  !! fit file missing or incomplete -- see /tmp/drfit_fit_${sample}_${wp}_${step}_${method}.log"
          continue
        fi

        # ---- Stage 3: plot + read-back verification ---------------------------------------------
        log "Stage 3 [${sample}/${wp}/step${step}/${method}]: plots + read-back check"
        root -l -b -q "plot_dr_correction_fits.cxx+(\"${sample}\", ${WPF}, ${step}, \"${method}\")" \
          > "/tmp/drfit_plot_${sample}_${wp}_${step}_${method}.log" 2>&1 || true

        NPNG=$(find "${MDIR}" -maxdepth 1 -name "step${step}_dr_fit_${method}_*.png" 2>/dev/null | wc -l)
        if [[ "${NPNG}" -lt 5 ]]; then
          ARTEFACT_FAILURES+=("only ${NPNG}/5 PNGs in ${MDIR}")
          log "  !! only ${NPNG} PNGs (expected 5: 4 pair-pT bins + inclusive)"
        fi
        # PERSISTENCE is a hard failure (the compiled-TF1 read-back trap); FLATNESS beyond Rp is
        # a property of the chosen function and is reported, not enforced -- `expo` is expected
        # to miss it, and that is exactly the information the method comparison needs.
        if [[ ! -f "${MDIR}readback_check.txt" ]]; then
          ARTEFACT_FAILURES+=("no readback_check.txt in ${MDIR}")
        else
          if ! grep -qE "^checked [0-9]+ functions, 0 FAILED persistence" "${MDIR}readback_check.txt"; then
            ARTEFACT_FAILURES+=("READ-BACK PERSISTENCE FAILURE in ${MDIR}readback_check.txt")
            log "  !! READ-BACK FAILURE -- persisted functions do not evaluate correctly"
          fi
          FLATLINE=$(grep "^flatness beyond Rp" "${MDIR}readback_check.txt" || true)
          [[ -n "${FLATLINE}" ]] && log "  ${FLATLINE}"
        fi

        # Summary line per method. The INCLUSIVE chi2/ndf and the median over sane-plateau cells
        # are the honest comparators; a plain mean is hostage to the 10k-overlay cells whose
        # plateau is not measurable at all (see fit_report.txt).
        if [[ -f "${MDIR}fit_report.txt" ]]; then
          INCL=$(grep "INCLUSIVE cell" "${MDIR}fit_report.txt" | tail -1 | awk -F'= ' '{print $NF}')
          SANE=$(grep "cells with |plateau-1|" "${MDIR}fit_report.txt" | tail -1 | sed 's/.*<= [0-9.]*: //')
          CHI2_SUMMARY+=("$(printf '%-10s %-7s step%-2s %-18s incl chi2/ndf = %-8s | sane cells: %s' \
                            "${sample}" "${wp}" "${step}" "${method}" "${INCL:-n/a}" "${SANE:-n/a}")")
        fi
      done
    done
  done
done

echo
echo "================================ SUMMARY ================================"
printf '%s\n' "${CHI2_SUMMARY[@]-}" | sed '/^$/d'
echo
if [[ ${#GUARD_FAILURES[@]} -gt 0 ]]; then
  echo "PLATEAU GUARD FAILED (FULL samples) for: ${GUARD_FAILURES[*]}"
  echo "  -> the fits above were produced with the explicit override; they are for inspection."
  echo "     See <plot base>/step<N>_dr_fit/plateau_guard_report.txt for the offending cells."
fi
if [[ ${#ARTEFACT_FAILURES[@]} -gt 0 ]]; then
  echo "ARTEFACT FAILURES:"
  printf '  %s\n' "${ARTEFACT_FAILURES[@]}"
fi
echo "========================================================================="

if [[ ${#ARTEFACT_FAILURES[@]} -gt 0 ]]; then exit 3; fi
if [[ ${#GUARD_FAILURES[@]}   -gt 0 ]]; then exit 2; fi
log "all stages complete, guard passed everywhere"
exit 0
