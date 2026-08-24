#!/usr/bin/env bash
# =============================================================================================
# MC TRIGGER-EFFICIENCY CLOSURE PIPELINE
# (docs/tracking/mc_trig_eff_closure.md; executes mc_trigger_efficiency.md Remaining Work 8)
#
# Closes the pp24 2mu4 trigger correction on the MC sample that carries BOTH an unbiased
# denominator and the per-pair trigger decision:
#
#     C = Sum_{2mu4} w_MC / [eps_MC(1) eps_MC(2) eps_dR(dR)]  /  Sum_{all} w_MC   = 1 ?
#
# The applied single-muon efficiency is the MC one (user, 2026-08-18; closure doc D4), so the test
# is SELF-CONTAINED and its residual is the dR correction alone.
#
# FOUR CELL-GROUPING APPROACHES (user, 2026-08-24;
# docs/tracking/mc_trigeff_dr_binning_approaches.md), each into its own output file and its own
# plot subdirectory. They differ ONLY in how the (pair pT, pair eta) plane is partitioned before
# the Step-3 fit:
#   nocorr                    8 pair-pT x 9 pair-eta = 72 cells   -- the un-merged reference
#   nocorr_ptmerge            7 x 9 = 63   -- last two pair-pT bins combined. 8-BIN AXIS ONLY.
#   nocorr_etamerge           8 x 3 = 24   -- pair eta combined into the 3 detector regions
#                                            (negative endcap / barrel / positive endcap)
#   nocorr_etamerge_ptmerge   7 x 3 = 21   -- both. 8-BIN AXIS ONLY.
# The DELIVERED correction is the SAME cascade in all four -- expo fit -> polyu fit -> `interp`
# interpolation -> the raw measured bins -- so what differs between the closures is the grouping.
# All three fit methods are therefore needed in EVERY mode. The un-merged reference additionally
# draws the two parametric forms as separate series, into a `separate_fit_forms/` subdirectory.
#
# EVERY closure is binned on the pp24 CROSS-SECTION's binning (ParamsSet::pT_bins_150, 15 log bins
# 8-150 GeV, x the 9 pair-eta panels), whatever its correction cells are -- doc D8. That is what
# makes the four comparable, and it is why Stage 4 can overlay them in one figure.
#
# STAGES
#   Stage 0  compile the two macros ONCE (a race on a shared _cxx.so is the real hazard)
#   Stage 1  UPSTREAM, 4-bin variant only (REGEN_4BIN=1). The canonical 8-bin Step-3 fits are
#            produced by the trigger-efficiency chain itself (run_mc_trigeff_round7.sh +
#            run_dr_correction_fits.sh) and are NOT rebuilt here. The 4-bin variant is:
#              (a) FillMCTrigEffHists Step 3   -> ..._pt4bin_step3.root
#              (b) plot_mc_trig_eff            -> dr_correction_plateaus_..._pt4bin.root
#              (c) fit_dr_corrections, OPPOSITE SIGN, plateau mode `nocorr`, per method
#            It has to be regenerated rather than reused: the `_pt4bin` outputs on disk predate
#            both the per-sign booking and the no-plateau-correction fit variant
#            (mc_trigger_efficiency.md R24b), so neither input this closure needs exists in them.
#   Stage 2  FillMCTrigEffClosure   -> mc_trig_eff_closure_<label><wp><ptbin><mode>.root
#   Stage 3  plot_mc_trig_eff_closure -> <plot base>/closure/<mode dir>/*.png
#            (2 PNGs per WP x binning x mode; the un-merged mode adds 2 more in
#             <mode dir>/separate_fit_forms/)
#   Stage 4  plot_mc_trig_eff_closure_compare -> <plot base>/closure/approach_comparison/*.png
#            the 4 approaches overlaid with the no-trigger series (5 lines). Needs ALL FOUR
#            approaches, so it runs for the canonical 8-bin binning only.
#
# WHY ARTEFACT VALIDATION AND NOT EXIT CODES: a ROOT macro that throws still exits 0 (the
# exception aborts the interpreter after ROOT has decided the batch job "ran"). Every stage is
# therefore checked by looking at what it actually produced.
#
# Usage:
#   ./run_mc_trigeff_closure.sh                      # both WPs, both binnings, no upstream regen
#   REGEN_4BIN=1 ./run_mc_trigeff_closure.sh         # ... and rebuild the 4-bin upstream first
#   WPS=tight PTBINS=8 ./run_mc_trigeff_closure.sh
#
# Env vars:
#   WPS="tight medium"    working points (registry: Analysis/docs/muon_wp_registry.md)
#   PTBINS="8 4"          pair-pT binnings; 4 sets MCTRIGEFF_PAIRPT_4BIN
#   MODES="nocorr nocorr_ptmerge nocorr_etamerge nocorr_etamerge_ptmerge"
#                         cell-grouping approaches. The two pair-pT-merging ones are skipped for
#                         PTBINS=4 (the merge is undefined there and THROWS).
#   SKIP_COMPARE=1        skip Stage 4 (the 4-approach overlay)
#   REGEN_4BIN=0|1        rebuild the 4-bin upstream (Stage 1). Costs ~4 min.
#   SKIP_FILL=1           reuse the closure ROOT files, re-make the plots only
# =============================================================================================
set -Eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
PLOT_DIR="${ANALYSIS_DIR}/plotting_codes/trig_effcy/mc_based"
LOG_DIR="${SCRIPT_DIR}/logs_closure"
mkdir -p "${LOG_DIR}"

SAMPLE="pp_full"          # pp only -- the 2mu4 PRODUCT weight; Pb+Pb is the mu4 union (doc Scope)
WPS="${WPS:-tight medium}"
PTBINS="${PTBINS:-8 4}"
REGEN_4BIN="${REGEN_4BIN:-0}"
SKIP_FILL="${SKIP_FILL:-0}"
# MUST match DrCorrCascadeMethods() in Utilities/DrCorrectionCascadeEvaluator.h -- the delivered
# cascade loads ALL THREE in EVERY mode (expo -> polyu_fixedRp -> interp), so all three fit files
# have to exist before any fill.
METHODS="expo polyu_fixedRp interp"
MODES="${MODES:-nocorr nocorr_ptmerge nocorr_etamerge nocorr_etamerge_ptmerge}"
SKIP_COMPARE="${SKIP_COMPARE:-0}"

MC_DIR="/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample"
PLOT_BASE="/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency"
LABEL="pp24_full"
# The DATA tag-and-probe turn-on directory -- MUST mirror kDataPPFitTmpl in FillMCTrigEffClosure.cxx
# (the C++ builds the real path; this shell only validates its freshness).
DATA_FIT_DIR="/usatlas/u/yuhanguo/usatlasdata/dimuon_data/pp_2024/trg_effcy_pT_fitting_to_erf_plus_log"

now() { date '+%F %T'; }
log() { echo "[$(now)] $*"; }
fail() { echo "[$(now)] ERROR: $*" >&2; exit 1; }
val()  { [[ -s "$1" ]] || fail "missing/empty artefact: $1"; }

# A compiled artefact is only trustworthy if the compile actually SUCCEEDED and the .so is newer
# than its source. `.L X.cxx+` leaves the PREVIOUS .so in place when the build fails, so checking
# only that the file exists validates last run's binary -- the documented "stale .so" anti-pattern.
val_compiled() {   # $1 = source, $2 = .so, $3 = compile log
  val "$2"
  # ACLiC's own banner is "Error in <ACLiC>: Compilation failed!" -- no colon straight after
  # "Error", so a pattern anchored on "error:" alone misses the very case this check exists for.
  grep -qE "(^| )[Ee]rror:|Error in <ACLiC>|Compilation failed|FATAL|fatal error" "$3" \
    && fail "compile errors in $3"
  # Compared against the HEADERS too: ACLiC leaves the previous .so in place when a build fails,
  # so after a header-only edit the stale .so is still newer than the .cxx and would pass.
  local src
  for src in "${@:4}" "$1"; do
    [[ -e "$src" ]] || continue
    [[ "$2" -nt "$src" ]] || fail "$2 is OLDER than $src -- the ACLiC build did not run (see $3)"
  done
  return 0
}

# An output built BEFORE one of its inputs is stale, and nothing downstream can tell. The concurrent
# trigger-efficiency session rewrites the dR fit files, so this is a live hazard, not a hypothetical.
val_fresh() {   # $1 = output, $2.. = inputs
  local out="$1"; shift
  local f
  for f in "$@"; do
    [[ -e "$f" ]] || fail "input does not exist: $f"
    [[ "$out" -nt "$f" ]] || fail "STALE OUTPUT: $out is older than its input $f
  -> the inputs were rewritten after it (a concurrent re-fit?). Re-run without SKIP_FILL."
  done
}

# `set -u` is deliberately NOT used: ~/setup.sh references unbound variables and dies under it.
log "sourcing ROOT environment"
source ~/setup.sh >/dev/null 2>&1 || fail "cannot source ~/setup.sh"
command -v root >/dev/null 2>&1 || fail "root not on PATH after sourcing ~/setup.sh"

wp_flag()   { [[ "$1" == "tight" ]] && echo "true" || echo "false"; }
wp_suffix() { [[ "$1" == "tight" ]] && echo ""     || echo "_medium_wp"; }
# Plot-tree suffix, mirroring DrCorrOutTag in dr_correction_sample_cfg.h. The C++ builds the real
# paths; this shell only VALIDATES them -- every time the two constructions have drifted apart the
# driver has reported real files as missing (see run_dr_correction_fits.sh).
plot_tag() { local t=""; [[ "$1" == "4" ]] && t="_pt4bin"; [[ "$2" == "tight" ]] || t="${t}_medium"; echo "$t"; }
ptbin_suf() { [[ "$1" == "4" ]] && echo "_pt4bin" || echo ""; }
# File token and plot subdirectory of an eps_dR variant. These MIRROR DrCorrPlateauModeTag /
# DrCorrPlateauModeDir in dr_correction_sample_cfg.h -- the C++ builds the real paths, this shell
# only validates them, and every time the two constructions drifted the driver reported real files
# as missing.
mode_tag() { case "$1" in
    nocorr)                  echo "_nocorr";;
    nocorr_ptmerge)          echo "_nocorr_ptmerge";;
    nocorr_etamerge)         echo "_nocorr_etamerge";;
    nocorr_etamerge_ptmerge) echo "_nocorr_etamerge_ptmerge";;
    *) fail "unknown eps_dR mode '$1'";; esac; }
mode_dir() { case "$1" in
    nocorr)                  echo "no_plateau_correction";;
    nocorr_ptmerge)          echo "no_plateau_correction_last2ptbins_merged";;
    nocorr_etamerge)         echo "no_plateau_correction_paireta_merged";;
    nocorr_etamerge_ptmerge) echo "no_plateau_correction_paireta_merged_last2ptbins_merged";;
    *) fail "unknown eps_dR mode '$1'";; esac; }
# Does this mode combine the last two pair-pT bins? Those are the 8-bin-axis-only ones.
mode_merges_pt() { case "$1" in nocorr_ptmerge|nocorr_etamerge_ptmerge) return 0;; *) return 1;; esac; }

# ---- Stage 0: compile ONCE ------------------------------------------------------------------
# `root -q 'X.cxx+'` with no argument compiles AND RUNS X() with its defaults; `-e '.L X.cxx+'`
# only builds.
log "Stage 0: compiling the macros with ACLiC (once)"
( cd "${RDF_DIR}"  && root -l -b -q -e '.L FillMCTrigEffClosure.cxx+' ) >"${LOG_DIR}/compile_fill.log" 2>&1
val_compiled "${RDF_DIR}/FillMCTrigEffClosure.cxx" "${RDF_DIR}/FillMCTrigEffClosure_cxx.so" \
             "${LOG_DIR}/compile_fill.log" \
             "${ANALYSIS_DIR}/Utilities/MCTrigEffPairSelection.h" \
             "${ANALYSIS_DIR}/Utilities/SingleMuEffEvaluator.h" \
             "${ANALYSIS_DIR}/Utilities/DrCorrectionCascadeEvaluator.h" \
             "${PLOT_DIR}/dr_correction_cell_groups.h" \
             "${ANALYSIS_DIR}/Utilities/MCTrigEffPairPtBinning.h" \
             "${ANALYSIS_DIR}/RDFBasedHistFilling/CommonEffcyConfig.h" \
             "${ANALYSIS_DIR}/MuonObjectsParamsAndHelpers/ParamsSet.h" \
             "${PLOT_DIR}/dr_correction_apply.h" \
             "${PLOT_DIR}/dr_correction_sample_cfg.h"
( cd "${PLOT_DIR}" && root -l -b -q -e '.L plot_mc_trig_eff_closure.cxx+' ) >"${LOG_DIR}/compile_plot.log" 2>&1
val_compiled "${PLOT_DIR}/plot_mc_trig_eff_closure.cxx" "${PLOT_DIR}/plot_mc_trig_eff_closure_cxx.so" \
             "${LOG_DIR}/compile_plot.log" \
             "${PLOT_DIR}/dr_correction_apply.h" "${PLOT_DIR}/dr_correction_ratio.h" \
             "${PLOT_DIR}/dr_correction_sample_cfg.h" \
             "${ANALYSIS_DIR}/Utilities/MCTrigEffPairPtBinning.h" \
             "${ANALYSIS_DIR}/Utilities/CommonLogYRange.h" \
             "${ANALYSIS_DIR}/RDFBasedHistFilling/CommonEffcyConfig.h"
( cd "${PLOT_DIR}" && root -l -b -q -e '.L plot_mc_trig_eff_closure_compare.cxx+' ) \
  >"${LOG_DIR}/compile_compare.log" 2>&1
val_compiled "${PLOT_DIR}/plot_mc_trig_eff_closure_compare.cxx" \
             "${PLOT_DIR}/plot_mc_trig_eff_closure_compare_cxx.so" \
             "${LOG_DIR}/compile_compare.log" \
             "${PLOT_DIR}/dr_correction_apply.h" "${PLOT_DIR}/dr_correction_ratio.h" \
             "${PLOT_DIR}/dr_correction_sample_cfg.h" \
             "${ANALYSIS_DIR}/Utilities/CommonLogYRange.h" \
             "${ANALYSIS_DIR}/RDFBasedHistFilling/CommonEffcyConfig.h"
if [[ "${REGEN_4BIN}" == "1" ]]; then
  ( cd "${RDF_DIR}"  && root -l -b -q -e '.L FillMCTrigEffHists.cxx+' ) >"${LOG_DIR}/compile_hists.log" 2>&1
  ( cd "${PLOT_DIR}" && root -l -b -q -e '.L plot_mc_trig_eff.cxx+' )   >"${LOG_DIR}/compile_trigeff.log" 2>&1
  ( cd "${PLOT_DIR}" && root -l -b -q -e '.L fit_dr_corrections.cxx+' ) >"${LOG_DIR}/compile_fit.log" 2>&1
fi
log "  compiled OK"

# ---- Stage 1: 4-bin upstream ----------------------------------------------------------------
if [[ "${REGEN_4BIN}" == "1" && " ${PTBINS} " == *" 4 "* ]]; then
  export MCTRIGEFF_PAIRPT_4BIN=1
  for wp in ${WPS}; do
    WPF="$(wp_flag "$wp")"; WPS_SUF="$(wp_suffix "$wp")"
    log "Stage 1a [4-bin/${wp}]: Step-3 fill"
    ( cd "${RDF_DIR}" && root -l -b -q "FillMCTrigEffHists.cxx+(\"${SAMPLE}\", true, ${WPF})" ) \
      >"${LOG_DIR}/regen4_fill_${wp}.log" 2>&1
    val "${MC_DIR}/mc_trig_eff_hists_${LABEL}${WPS_SUF}_pt4bin_step3.root"
    log "Stage 1b [4-bin/${wp}]: plateau measurement"
    ( cd "${PLOT_DIR}" && root -l -b -q "plot_mc_trig_eff.cxx+(\"${SAMPLE}\", ${WPF})" ) \
      >"${LOG_DIR}/regen4_plateau_${wp}.log" 2>&1
    val "${MC_DIR}/dr_correction_plateaus_${LABEL}${WPS_SUF}_pt4bin.root"
    for mode in ${MODES}; do
      # The pair-pT merge is defined on the canonical 8-bin axis only -- the 4-bin variant already
      # combines the top cells by construction and the C++ THROWS if asked.
      if mode_merges_pt "$mode"; then continue; fi
      MTAG="$(mode_tag "$mode")"
      for m in ${METHODS}; do
        log "Stage 1c [4-bin/${wp}/${mode}/${m}]: opposite-sign fit, no plateau correction"
        ( cd "${PLOT_DIR}" && root -l -b -q \
            "fit_dr_corrections.cxx+(\"${SAMPLE}\", ${WPF}, 3, \"${m}\", false, \"os\", \"${mode}\")" ) \
          >"${LOG_DIR}/regen4_fit_${wp}_${mode}_${m}.log" 2>&1
        val "${MC_DIR}/dr_correction_fits_${LABEL}${WPS_SUF}_pt4bin_step3_${m}_os${MTAG}.root"
      done
    done
  done
  unset MCTRIGEFF_PAIRPT_4BIN
fi

# ---- Stages 2 + 3 ----------------------------------------------------------------------------
for pb in ${PTBINS}; do
  if [[ "$pb" == "4" ]]; then export MCTRIGEFF_PAIRPT_4BIN=1; else unset MCTRIGEFF_PAIRPT_4BIN; fi
  PBSUF="$(ptbin_suf "$pb")"
  for wp in ${WPS}; do
    WPF="$(wp_flag "$wp")"; WPS_SUF="$(wp_suffix "$wp")"
    for mode in ${MODES}; do
      # The pair-pT merge is defined on the canonical 8-bin axis only
      # (dr_correction_cell_groups.h throws otherwise), so skip it loudly rather than letting the
      # macro die inside ROOT. The pair-eta merge is orthogonal and runs on either axis.
      if [[ "$pb" == "4" ]] && mode_merges_pt "$mode"; then
        log "Stage 2/3 [4-bin/${wp}/${mode}]: SKIPPED -- combining the last two pair-pT bins is "\
"defined for the canonical 8-bin pair-pT axis only"
        continue
      fi
      MTAG="$(mode_tag "$mode")"; MDIR="$(mode_dir "$mode")"
      OUT="${MC_DIR}/mc_trig_eff_closure_${LABEL}${WPS_SUF}${PBSUF}${MTAG}.root"

      # Every fit file the fill will read must exist BEFORE the fill, or the macro throws inside
      # ROOT and (exit code 0) looks like success until the artefact check three lines later.
      # ALL THREE methods, in EVERY mode: the delivered cascade loads expo, polyu_fixedRp and
      # interp and routes each cell to the first one that is accepted.
      for m in ${METHODS}; do
        f="${MC_DIR}/dr_correction_fits_${LABEL}${WPS_SUF}${PBSUF}_step3_${m}_os${MTAG}.root"
        [[ -s "$f" ]] || fail "no opposite-sign '${mode}' fit for ${pb}-bin/${wp}/${m}:
    ${f}
  Produce it with run_dr_correction_fits.sh (SIGNS=os PLATEAU_MODES=${mode}) for the 8-bin
  binning, or re-run this script with REGEN_4BIN=1 for the 4-bin variant."
      done

      if [[ "${SKIP_FILL}" == "1" ]]; then
        log "Stage 2 [${pb}-bin/${wp}/${mode}]: SKIPPED (SKIP_FILL=1)"
      else
        log "Stage 2 [${pb}-bin/${wp}/${mode}]: closure fill"
        rm -f "${OUT}"          # so the check below cannot pass on a stale file
        ( cd "${RDF_DIR}" && root -l -b -q \
            "FillMCTrigEffClosure.cxx+(\"${SAMPLE}\", ${WPF}, \"${mode}\")" ) \
          >"${LOG_DIR}/fill_${pb}bin_${wp}_${mode}.log" 2>&1
      fi
      val "${OUT}"
      # FRESHNESS, not just existence: the closure must be newer than every efficiency input it was
      # weighted by. Its own provenance TNamed records those files' mtimes for the same reason.
      # EVERY input that enters the weight, not just the dR fits: the single-muon turn-ons (the
      # APPLIED eps_MC and the eps^nc_data of the diagnostic numerator) are rewritten by the
      # concurrent trigger-efficiency session too, and they were outside the original gate.
      FRESH_INPUTS=("${MC_DIR}/mc_trig_eff_hists_${LABEL}${WPS_SUF}${PBSUF}_step3.root"
                    "${MC_DIR}/muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_full.root"
                    "${DATA_FIT_DIR}/single_mu_effcy_pT_fit${WPS_SUF}.root"
                    "${MC_DIR}/single_mu_effcy_pT_fit_mc${WPS_SUF}.root")
      for m in ${METHODS}; do
        FRESH_INPUTS+=("${MC_DIR}/dr_correction_fits_${LABEL}${WPS_SUF}${PBSUF}_step3_${m}_os${MTAG}.root")
      done
      val_fresh "${OUT}" "${FRESH_INPUTS[@]}"
      # The raw-bin fallback is a TEMPORARY PLACEHOLDER (doc §3.3): surface it every run, so it
      # cannot quietly become permanent.
      grep -h "RAW-BIN PLACEHOLDER\|on the RAW measured bins" \
        "${LOG_DIR}/fill_${pb}bin_${wp}_${mode}.log" 2>/dev/null | sed 's/^/  /'

      log "Stage 3 [${pb}-bin/${wp}/${mode}]: closure plots"
      PDIR="${PLOT_BASE}/mc_based$(plot_tag "$pb" "$wp")/closure/${MDIR}"
      PNGS=(closure_pair_pt_all_opposite_sign.png closure_pair_pt_single_b_signal_cuts.png)
      # The un-merged reference produces a SECOND figure set -- the two parametric forms as
      # separate series -- in its own subdirectory (user, 2026-08-24). Mirrors kPerFormSubdir in
      # plot_mc_trig_eff_closure.cxx; the C++ builds the real path, this only validates it.
      PDIRS=("${PDIR}")
      [[ "$mode" == "nocorr" ]] && PDIRS+=("${PDIR}/separate_fit_forms")
      # Removed BEFORE the macro runs, for the same reason ${OUT} is: otherwise a plot pass that
      # dies inside ROOT (which still exits 0) validates on the previous run's PNGs.
      for d in "${PDIRS[@]}"; do for f in "${PNGS[@]}"; do rm -f "${d}/${f}"; done; done
      ( cd "${PLOT_DIR}" && root -l -b -q \
          "plot_mc_trig_eff_closure.cxx+(\"${SAMPLE}\", ${WPF}, \"${mode}\")" ) \
        >"${LOG_DIR}/plot_${pb}bin_${wp}_${mode}.log" 2>&1
      for d in "${PDIRS[@]}"; do for f in "${PNGS[@]}"; do
        val "${d}/${f}"
        val_fresh "${d}/${f}" "${OUT}"
      done; done
      log "  $(( 2 * ${#PDIRS[@]} )) PNGs in ${PDIR}"
      grep -h "inclusive closure" "${LOG_DIR}/fill_${pb}bin_${wp}_${mode}.log" 2>/dev/null | sed 's/^/  /'
    done
  done
done
unset MCTRIGEFF_PAIRPT_4BIN

# ---- Stage 4: the 4-approach overlay ---------------------------------------------------------
# Needs ALL FOUR approaches, two of which are defined for the canonical 8-bin pair-pT axis only,
# so it runs for the 8-bin binning alone. Its per-approach points must EQUAL the ones in each
# approach's own subdirectory -- same fills, same cascade, only overlaid -- and the macro itself
# enforces the precondition that makes the overlay legitimate: the four no-trigger denominators
# are compared bin by bin and it THROWS if they differ.
NEEDED_MODES="nocorr nocorr_ptmerge nocorr_etamerge nocorr_etamerge_ptmerge"
HAVE_ALL=1
for need in ${NEEDED_MODES}; do
  [[ " ${MODES} " == *" ${need} "* ]] || HAVE_ALL=0
done
if [[ "${SKIP_COMPARE}" == "1" ]]; then
  log "Stage 4: SKIPPED (SKIP_COMPARE=1)"
elif [[ "${HAVE_ALL}" != "1" || " ${PTBINS} " != *" 8 "* ]]; then
  log "Stage 4: SKIPPED -- the overlay needs all four approaches on the 8-bin axis (MODES='${MODES}', PTBINS='${PTBINS}')"
else
  CMP_PNGS=(closure_compare_pair_pt_all_opposite_sign.png
            closure_compare_pair_pt_single_b_signal_cuts.png)
  for wp in ${WPS}; do
    WPF="$(wp_flag "$wp")"; WPS_SUF="$(wp_suffix "$wp")"
    CDIR="${PLOT_BASE}/mc_based$(plot_tag "8" "$wp")/closure/approach_comparison"
    log "Stage 4 [${wp}]: the 4-approach overlay"
    for f in "${CMP_PNGS[@]}"; do rm -f "${CDIR}/${f}"; done
    ( cd "${PLOT_DIR}" && root -l -b -q \
        "plot_mc_trig_eff_closure_compare.cxx+(\"${SAMPLE}\", ${WPF})" ) \
      >"${LOG_DIR}/compare_${wp}.log" 2>&1
    CMP_INPUTS=()
    for mode in ${NEEDED_MODES}; do
      CMP_INPUTS+=("${MC_DIR}/mc_trig_eff_closure_${LABEL}${WPS_SUF}$(mode_tag "$mode").root")
    done
    for f in "${CMP_PNGS[@]}"; do
      val "${CDIR}/${f}"
      val_fresh "${CDIR}/${f}" "${CMP_INPUTS[@]}"
    done
    log "  2 PNGs in ${CDIR}"
    grep -h "inclusive closure" "${LOG_DIR}/compare_${wp}.log" 2>/dev/null | sed 's/^/  /'
  done
fi

log "════ MC trigger-efficiency closure DONE ════"
