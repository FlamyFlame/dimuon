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
#   Stage 2  GUARD + FIT         fit_dr_corrections(sample, wp, step, method, sign, plateau mode)
#              reads that ROOT file (never a .txt / .md), enforces |plateau-1| <= 0.15 on a FULL
#              production and divides each cell's curve by ITS OWN plateau -- or, in the
#              no-plateau-correction mode, applies neither and fits a free baseline instead
#   Stage 3  PLOT + READ-BACK    plot_dr_correction_fits(sample, wp, step, method, mode, pmode)
#              re-opens the fit file(s) in a SEPARATE process, draws measured vs fitted -- 1 PNG
#              per pair-pT bin, 1 subplot per pair-eta bin -- and verifies every persisted
#              function outside its fit range
#
# PLATEAU MODE. The dR fit runs in one or both of two modes, and the mode is the TOP level of
# the plot tree so the two can never interleave:
#   step<N>_dr_fit/plateau_corrected/<method>/<sign mode>/       NOMINAL. Each cell's curve is
#              divided by ITS OWN large-dR plateau and the fitted shape tends to 1.
#   step<N>_dr_fit/no_plateau_correction/<method>/<sign mode>/   The RAW efficiency is fitted with
#              a FREE baseline C, so the baseline comes from the dR < 1 data itself and the
#              [2, 3.5] plateau window never enters the fit. Motivation (user, 2026-08-11): the
#              inverse-weighted dR distribution carries structure out to large dR, worst in the
#              pair-eta bins enclosing the detector gap, so a far-away plateau may be the wrong
#              baseline for the small-dR region the correction is about.
#   step<N>_dr_fit/no_plateau_correction_last2ptbins_merged/<method>/<sign mode>/
#              The SAME raw fit with the LAST TWO pair-pT bins MERGED into one cell (user,
#              2026-08-17), so this tree has ONE canvas fewer and its top one covers
#              p_T^pair in [74.24, 150) GeV. The top two cells of the 8-bin log axis run past where
#              the sample has yield. It is NOT a new binning -- the two bins are projected together
#              at the fit stage; see dr_correction_cell_groups.h. 8-BIN AXIS ONLY: it is skipped,
#              with a printed note, when MCTRIGEFF_PAIRPT_4BIN is set. With opposite-sign pairs and
#              the `expo` form this is the variant the pp24 crossx application uses.
#   step<N>_dr_fit/no_plateau_correction_paireta_merged/<method>/<sign mode>/
#              The SAME raw fit with the pair-eta bins MERGED into THREE SIGN-INDEPENDENT
#              |eta^pair| BINS -- |eta| < 1.0 (barrel), 1.0 <= |eta| < 2.0, 2.0 <= |eta| < the source axis's top edge (2.2 since 2026-09-07,
#                     tracking ParamsSet::pair_eta_fiducial_max)
#              (user, 2026-08-24; SUPERSEDED 2026-09-03 -- the dR correlation barely depends on
#              the SIGN of pair eta, but the |eta|>2-vs-<2 split inside the endcap is much bigger
#              than any negative/positive asymmetry) -- so every canvas here carries 3 panels
#              instead of 9. The 9-bin pair-eta grid is the CROSS-SECTION's presentation binning,
#              never chosen for eps_dR's statistics; the three bins triple the pairs per fit while
#              keeping the one physically motivated distinction. The barrel group is one
#              contiguous source range; each forward group FOLDS its negative- and positive-eta
#              source bins together. Also NOT a new binning -- the source bins are projected
#              together at the fit stage and every |eta| boundary is looked up as a symmetric
#              pair of existing edges; see dr_correction_cell_groups.h. It is ORTHOGONAL to the
#              pair-pT axis, so it runs on the 4-bin variant too.
#   step<N>_dr_fit/no_plateau_correction_paireta_merged_last2ptbins_merged/<method>/<sign mode>/
#              Both merges at once. Because it merges pair pT as well, the 8-BIN-AXIS-ONLY
#              restriction above applies to it too and it is skipped under MCTRIGEFF_PAIRPT_4BIN.
# STEP 4 RUNS THE NOMINAL MODE ONLY (user: "ignore step4, focus on step3 for now"), but its output
# still lands under plateau_corrected/ so both step trees have the same shape.
#
# SIGN SERIES. One run produces up to three fitted series per (sample, WP, step, method):
# sign-integrated (the NOMINAL correction), same sign and opposite sign. They differ only in the
# input histograms and plateaus; the guard is FATAL only for the sign-integrated series of a FULL
# production (each sign holds a fraction of the pairs, so per-sign failures are a statistics
# statement -- see fit_dr_corrections.cxx). The plots come out as two sets:
#   step<N>_dr_fit/<plateau mode>/<method>/sign_intgr/  the sign-integrated series alone
#   step<N>_dr_fit/<plateau mode>/<method>/sign_sepr/   same sign and opposite sign overlaid
# A sample without per-sign inputs (e.g. an overlay filled before the per-sign booking) skips its
# sign-separated series with a printed note and still produces the sign-integrated one.
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
#   PLATEAU_MODES="corr nocorr nocorr_ptmerge nocorr_etamerge nocorr_etamerge_ptmerge"
#                                    plateau modes to run. `corr` = the nominal, plateau-
#                                    normalized fit; `nocorr` = the raw fit with a free baseline;
#                                    `nocorr_ptmerge` = the same with the last two pair-pT bins
#                                    merged (8-bin axis only); `nocorr_etamerge` = the same with
#                                    the pair-eta bins folded into 3 sign-independent |eta| bins
#                                    (no axis restriction); `nocorr_etamerge_ptmerge` = both merges
#                                    (8-bin axis only, because it merges pair pT). Step 4 is
#                                    restricted to `corr` whatever this says.
#   SIGNS="intgr ss os"              sign series to fit. `intgr` = sign-integrated (nominal),
#                                    `ss` = same sign, `os` = opposite sign. The sign_sepr plot
#                                    set needs BOTH ss and os.
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
# Default methods (user 2026-08-04): NOMINAL `expo`, BACKUP `polyu_fixedRp`, cross-check
# `interp`. The two power laws are DROPPED -- see MakeMethodCfg in fit_dr_corrections.cxx: the
# fitted exponent n rails to its lower limit 0.2 in many cells and n < 1 gives an INFINITE slope
# at Rp (a visible cusp), so they are not smooth. They remain constructible for reproducing old
# outputs, but nothing produces them by default.
METHODS="${METHODS:-expo polyu_fixedRp interp}"
# Sign series. "intgr" is the token for the sign-INTEGRATED (nominal) series; the C++ takes "" for
# it, so it is translated in sign_arg() below and never typed as an empty word here.
SIGNS="${SIGNS:-intgr ss os}"
# Plateau modes. The C++ takes "corr" / "nocorr"; the directory names and the fit-file token are
# built by dr_correction_sample_cfg.h (DrCorrPlateauModeDir / DrCorrPlateauModeTag) and only
# MIRRORED here, in pmode_dir/pmode_fsuf below, for the artefact checks.
PLATEAU_MODES="${PLATEAU_MODES:-corr nocorr nocorr_ptmerge nocorr_etamerge nocorr_etamerge_ptmerge}"
# The pair-pT-binning token must MIRROR Utilities/MCTrigEffPairPtBinning.h: the C++ writes
# ..._pt4bin... when MCTRIGEFF_PAIRPT_4BIN is set, and these artefact checks look the files up by
# name. When they disagreed, a perfectly good 4-bin run was reported as 18 "missing fit files".
PTBIN_SUF=""; [[ -n "${MCTRIGEFF_PAIRPT_4BIN:-}" ]] && PTBIN_SUF="_pt4bin"
PTBIN_DIR=""   # retired: the pair-pT variant is in the top-level base now
# Pair-pT binning token -- MUST match Utilities/MCTrigEffPairPtBinning.h::FileSuffix().
# The C++ builds the real filenames; this shell only VALIDATES them, and a name built in two
# places is exactly how the 4-bin pass came to report every artefact as missing while the
# files were in fact written correctly under their token.
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
#
# NO `-q` HERE, and that is not a style choice: `root -l -b -q` with no macro argument QUITS
# BEFORE READING STDIN, so the heredoc below was never executed and this function returned 0
# unconditionally -- i.e. every artefact check it backs silently PASSED, missing objects and all.
# Verified 2026-08-11: `root -l -b -q` fed `gSystem->Exit(7);` on stdin exits 0; `root -l -b`
# exits 7. Fixed 2026-08-11 together with the sign-series work.
root_file_has_objects() {
  local f="$1"; shift
  [[ -s "$f" ]] || return 1
  local names=("$@") checks=""
  for n in "${names[@]}"; do
    checks+="if (!fin->Get(\"${n}\")) { fin->Close(); gSystem->Exit(4); }"$'\n'
  done
  root -l -b <<EOF >/dev/null 2>&1
TFile *fin = TFile::Open("$f", "READ");
if (!fin || fin->IsZombie()) { gSystem->Exit(2); }
${checks}
fin->Close();
gSystem->Exit(0);
EOF
}

# sign token -> the C++ argument / the file-name suffix / the human wording used in report names.
# Mirrors dr_correction_sample_cfg.h (DrCorrSignText / DrCorrSignFileTag): the C++ builds the real
# names, this shell only VALIDATES them.
sign_arg()  { [[ "$1" == "intgr" ]] && echo "" || echo "$1"; }
sign_fsuf() { [[ "$1" == "intgr" ]] && echo "" || echo "_$1"; }
sign_rtag() {
  case "$1" in
    intgr) echo "" ;;
    ss)    echo "_same_sign" ;;
    os)    echo "_opposite_sign" ;;
    *) fail "unknown sign '$1' (use intgr | ss | os)" ;;
  esac
}
# Plateau-mode -> plot subdirectory / fit-file token. MUST mirror dr_correction_sample_cfg.h:
# the C++ builds the real names, this shell only VALIDATES what it wrote.
PMODE_TOKENS="corr | nocorr | nocorr_ptmerge | nocorr_etamerge | nocorr_etamerge_ptmerge"
pmode_dir() {
  case "$1" in
    corr)                    echo "plateau_corrected/" ;;
    nocorr)                  echo "no_plateau_correction/" ;;
    nocorr_ptmerge)          echo "no_plateau_correction_last2ptbins_merged/" ;;
    nocorr_etamerge)         echo "no_plateau_correction_paireta_merged/" ;;
    nocorr_etamerge_ptmerge) echo "no_plateau_correction_paireta_merged_last2ptbins_merged/" ;;
    *) fail "unknown plateau mode '$1' (use ${PMODE_TOKENS})" ;;
  esac
}
pmode_fsuf() {
  case "$1" in
    corr)                    echo "" ;;
    nocorr)                  echo "_nocorr" ;;
    nocorr_ptmerge)          echo "_nocorr_ptmerge" ;;
    nocorr_etamerge)         echo "_nocorr_etamerge" ;;
    nocorr_etamerge_ptmerge) echo "_nocorr_etamerge_ptmerge" ;;
    *) fail "unknown plateau mode '$1' (use ${PMODE_TOKENS})" ;;
  esac
}
pmode_text() {
  case "$1" in
    corr)                    echo "plateau-corrected" ;;
    nocorr)                  echo "no plateau correction" ;;
    nocorr_ptmerge)          echo "no plateau corr., last 2 pT bins merged" ;;
    nocorr_etamerge)         echo "no plateau corr., pair-eta folded into 3 |eta| bins" ;;
    nocorr_etamerge_ptmerge) echo "no plateau corr., pair-eta folded into 3 |eta| bins + last 2 pT bins merged" ;;
  esac
}

sign_text() {
  case "$1" in
    intgr) echo "sign-integrated" ;;
    ss)    echo "same sign" ;;
    os)    echo "opposite sign" ;;
  esac
}

# Number of pair-pT bins, READ FROM THE FIT FILE ITSELF (the x axis of h_step<N>_plateau) rather
# than hard-coded. The expected PNG count is that number + 1 (one canvas per pair-pT bin plus the
# inclusive one), so it follows the binning automatically: the literal 5 this check used to carry
# was the 4-bin variant's count and silently under-checked the 8-bin nominal (which makes 9).
fit_file_npt() {
  local f="$1" step="$2"
  [[ -s "$f" ]] || { echo 0; return; }
  root -l -b <<EOF 2>/dev/null | tail -1
TFile *fin = TFile::Open("$f", "READ");
TH2 *h = fin ? (TH2*)fin->Get("h_step${step}_plateau") : nullptr;
printf("%d\\n", h ? h->GetNbinsX() : 0);
gSystem->Exit(0);
EOF
}

wp_flag()   { [[ "$1" == "tight" ]] && echo "true" || echo "false"; }
wp_suffix() { [[ "$1" == "tight" ]] && echo ""     || echo "_medium_wp"; }
# retired: the WP is in the top-level base now, not a per-directory subdir
wp_dir()    { echo ""; }

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
# MUST mirror dr_correction_sample_cfg.h::DrCorrOutTag + the out_base strings: each variant is
# its OWN top-level tree, mc_based{_pt4bin}{_medium}. This shell only VALIDATES what the C++
# wrote, and every time the two constructions have drifted apart the driver has either reported
# real files as missing or (worse) deleted the other variant's outputs. Derive, never retype.
sample_plot_base() {
  local sample="$1" wp="$2" tag="${PTBIN_SUF}"
  [[ "$wp" != "tight" ]] && tag="${tag}_medium"
  case "$sample" in
    pp|pp_full) echo "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/mc_based${tag}/" ;;
    overlay)    echo "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pbpb_trigger_efficiency/mc_based${tag}/" ;;
    noovl)      echo "/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/r17663_no_overlay_trigger_efficiency/mc_based${tag}/" ;;
    *) fail "unknown sample '$sample'" ;;
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

  for wp in ${WPS}; do
    WPF="$(wp_flag "${wp}")"
    WPS_SUF="$(wp_suffix "${wp}")"
    WPD="$(wp_dir "${wp}")"
    PLOTBASE="$(sample_plot_base "${sample}" "${wp}")"
    PLATEAU_FILE="${MCDIR}dr_correction_plateaus_${LABEL}${WPS_SUF}${PTBIN_SUF}.root"

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

      # STEP 4 IS NOMINAL-MODE ONLY (user: ignore step 4 for now). Its plots still land under
      # plateau_corrected/ so the step-3 and step-4 trees have the same shape.
      STEP_PMODES="${PLATEAU_MODES}"
      if [[ "${step}" == "4" ]]; then
        if [[ " ${PLATEAU_MODES} " == *" corr "* ]]; then
          STEP_PMODES="corr"
        else
          log "step 4 runs the plateau-corrected mode only -- not in PLATEAU_MODES, skipping it"
          continue
        fi
      fi

      for pmode in ${STEP_PMODES}; do
        # THE PAIR-pT MERGE IS DEFINED FOR THE 8-BIN NOMINAL AXIS ONLY (the C++ throws on the
        # other one; dr_correction_cell_groups.h). Skipping it here with a note keeps a 4-bin
        # production running instead of filling the summary with expected failures. EVERY mode
        # that merges pair pT is covered, `nocorr_etamerge_ptmerge` included -- it merges the top
        # two pair-pT bins exactly as `nocorr_ptmerge` does, so the same restriction applies.
        # `nocorr_etamerge` ALONE is NOT skipped: the pair-eta merge is orthogonal to the pair-pT
        # axis and is defined on both variants.
        if [[ ( "${pmode}" == "nocorr_ptmerge" || "${pmode}" == "nocorr_etamerge_ptmerge" ) \
              && -n "${MCTRIGEFF_PAIRPT_4BIN:-}" ]]; then
          log "  skipping plateau mode ${pmode}: it merges the last two pair-pT bins, which is"
          log "        defined for the 8-bin pair-pT axis only, and MCTRIGEFF_PAIRPT_4BIN is set"
          log "        (that variant already merges the top cells)"
          continue
        fi
        PMDIR="$(pmode_dir  "${pmode}")"
        PMSUF="$(pmode_fsuf "${pmode}")"
        PMTEXT="$(pmode_text "${pmode}")"

      for method in ${METHODS}; do
        MDIR="${PLOTBASE}step${step}_dr_fit/${PMDIR}${PTBIN_DIR}${method}/${WPD}"

        # ---- Stage 2: guard + fit, once per SIGN SERIES ----------------------------------------
        FITTED_SIGNS=""      # signs whose fit file came out complete -> drive Stage 3
        for sgn in ${SIGNS}; do
          SARG="$(sign_arg   "${sgn}")"
          SFSUF="$(sign_fsuf "${sgn}")"
          SRTAG="$(sign_rtag "${sgn}")"
          STEXT="$(sign_text "${sgn}")"
          FIT_FILE="${MCDIR}dr_correction_fits_${LABEL}${WPS_SUF}${PTBIN_SUF}_step${step}_${method}${SFSUF}${PMSUF}.root"
          GUARD_REPORT="${PLOTBASE}step${step}_dr_fit/${PMDIR}${PTBIN_DIR}${WPD}plateau_guard_report${SRTAG}.txt"
          FITLOG="/tmp/drfit_fit_${sample}_${wp}_${step}_${method}_${sgn}_${pmode}.log"

          if [[ "${SKIP_FIT}" == "1" ]]; then
            log "Stage 2 [${sample}/${wp}/step${step}/${method}/${STEXT}/${PMTEXT}]: SKIPPED (SKIP_FIT=1)"
          else
            # deleted first so the artefact check below cannot pass on a stale file
            rm -f "${FIT_FILE}"
            log "Stage 2 [${sample}/${wp}/step${step}/${method}/${STEXT}/${PMTEXT}]: guard + fit"
            root -l -b -q "fit_dr_corrections.cxx+(\"${sample}\", ${WPF}, ${step}, \"${method}\", false, \"${SARG}\", \"${pmode}\")" \
              > "${FITLOG}" 2>&1 || true
          fi

          # The guard writes its verdict BEFORE it throws, so the verdict is readable either way.
          # Only the SIGN-INTEGRATED series has a fatal tier: that series is the nominal
          # deliverable, while one charge combination carries a fraction of the statistics, so a
          # per-sign violation is a statistics statement and is reported, not enforced
          # (fit_dr_corrections.cxx writes "reported, not enforced" as its verdict there).
          # The fatal tier exists only where the plateau is actually APPLIED. In the
          # no-plateau-correction mode nothing is normalized by it, so it cannot disqualify
          # anything and its report's verdict says "not applicable".
          GUARD_FAILED=0
          if [[ "${sgn}" == "intgr" && "${pmode}" == "corr" && -f "${GUARD_REPORT}" ]] \
             && grep -q "^verdict: FAIL" "${GUARD_REPORT}"; then
            GUARD_FAILED=1
          fi

          if [[ "${GUARD_FAILED}" == "1" ]]; then
            if [[ ! " ${GUARD_FAILURES[*]-} " == *" ${sample}/${wp}/step${step} "* ]]; then
              GUARD_FAILURES+=("${sample}/${wp}/step${step}")
            fi
            echo "==============================================================================="
            echo " PLATEAU GUARD FAILED: ${sample} (FULL production) / ${wp} / step ${step}"
            sed -n '/^FAILING cells/,$p' "${GUARD_REPORT}"
            echo "==============================================================================="
            if [[ "${STRICT_GUARD}" == "1" ]]; then
              log "  STRICT_GUARD=1 -> no fits produced for ${sample}/${wp}/step${step}"
              continue
            fi
            if [[ "${SKIP_FIT}" != "1" ]]; then
              log "  re-running WITH the explicit override so the plots exist for human inspection"
              root -l -b -q "fit_dr_corrections.cxx+(\"${sample}\", ${WPF}, ${step}, \"${method}\", true, \"${SARG}\", \"${pmode}\")" \
                > "${FITLOG}" 2>&1 || true
            fi
          fi

          if ! root_file_has_objects "${FIT_FILE}" "h_step${step}_chi2ndf" \
                                                   "h_step${step}_plateau"; then
            # A sign-separated series whose per-sign inputs do not exist yet is SKIPPED on purpose
            # by the macro (it says so on stdout) -- that is not an artefact failure, it is a
            # sample that has not been re-filled with the per-sign booking.
            if [[ "${sgn}" != "intgr" ]] && grep -q "SKIPPED:" "${FITLOG}" 2>/dev/null; then
              log "  note: no ${STEXT} inputs for ${sample} -- that series is skipped (see ${FITLOG})"
            else
              ARTEFACT_FAILURES+=("fit file missing/incomplete: ${FIT_FILE}")
              log "  !! fit file missing or incomplete -- see ${FITLOG}"
            fi
            continue
          fi
          FITTED_SIGNS="${FITTED_SIGNS} ${sgn}"

          # Summary line per method and series. The INCLUSIVE chi2/ndf and the median over
          # sane-plateau cells are the honest comparators; a plain mean is hostage to the
          # 10k-overlay cells whose plateau is not measurable at all (see fit_report.txt).
          if [[ -f "${MDIR}fit_report${SRTAG}.txt" ]]; then
            INCL=$(grep "INCLUSIVE cell" "${MDIR}fit_report${SRTAG}.txt" | tail -1 | awk -F'= ' '{print $NF}')
            # The "sane plateau" subset only exists in the nominal mode (a curve divided by a bad
            # plateau has an inflated chi2). Nothing is divided in the other mode, so the honest
            # comparator there is the statistic over ALL converged cells.
            if [[ "${pmode}" == "corr" ]]; then
              SANE=$(grep "cells with |plateau-1|" "${MDIR}fit_report${SRTAG}.txt" | tail -1 | sed 's/.*<= [0-9.]*: //')
            else
              SANE=$(grep "over all converged cells" "${MDIR}fit_report${SRTAG}.txt" | tail -1 | sed 's/.*converged cells: *//')
            fi
            # ... and the label must say WHICH subset the number is, or the two modes' summary
            # lines look comparable when they are not.
            SUBSET="sane-plateau cells"; [[ "${pmode}" == "corr" ]] || SUBSET="all fitted cells"
            CHI2_SUMMARY+=("$(printf '%-10s %-7s step%-2s %-22s %-16s %-16s incl chi2/ndf = %-8s | %s: %s' \
                              "${sample}" "${wp}" "${step}" "${PMTEXT}" "${method}" "${STEXT}" "${INCL:-n/a}" "${SUBSET}" "${SANE:-n/a}")")
          fi
        done

        # ---- Stage 3: plot + read-back verification, one directory per MODE ---------------------
        # sign_intgr needs the sign-integrated fits; sign_sepr overlays both charge combinations
        # and therefore needs BOTH per-sign fit files.
        MODES=""
        [[ " ${FITTED_SIGNS} " == *" intgr "* ]] && MODES="sign_intgr"
        if [[ " ${FITTED_SIGNS} " == *" ss "* && " ${FITTED_SIGNS} " == *" os "* ]]; then
          MODES="${MODES} sign_sepr"
        elif [[ " ${SIGNS} " == *" ss "* || " ${SIGNS} " == *" os "* ]]; then
          log "  note: no sign_sepr plots for ${sample}/${wp}/step${step}/${method} -- that needs"
          log "        BOTH the same-sign and the opposite-sign fits (have:${FITTED_SIGNS:-none})"
        fi

        for mode in ${MODES}; do
        # dR VIEWS (user, 2026-08-12). STEP 3 carries two views of the same fitted cells: the
        # DEFAULT dR in [0,1] -- the fit domain, the window the correction is applied in, written
        # at the top of <sign mode>/ -- and the REFERENCE dR in [0,2] one level down in dR0_2/,
        # which keeps the large-dR points that DEFINED the plateau visible. Step 4 was NOT
        # restructured (the user asked for step3_dr_fit only), so it runs the [0,2] view alone and
        # keeps its flat layout; the C++ suppresses the dR0_2/ level for step != 3.
        VIEWS="dr0_1 dr0_2"; [[ "${step}" == "3" ]] || VIEWS="dr0_2"
        for view in ${VIEWS}; do
          VIEWDIR=""; [[ "${step}" == "3" && "${view}" == "dr0_2" ]] && VIEWDIR="dR0_2/"
          log "Stage 3 [${sample}/${wp}/step${step}/${method}/${mode}/${PMTEXT}/${view}]: plots + read-back check"
          root -l -b -q "plot_dr_correction_fits.cxx+(\"${sample}\", ${WPF}, ${step}, \"${method}\", \"${mode}\", \"${pmode}\", \"${view}\")" \
            > "/tmp/drfit_plot_${sample}_${wp}_${step}_${method}_${mode}_${pmode}_${view}.log" 2>&1 || true

          # Expected PNG count DERIVED from the binning in the fit file: one canvas per pair-pT
          # bin plus the inclusive one. Never a literal -- the literal 5 that used to be here was
          # the 4-bin variant's count and passed silently on the 8-bin nominal, which makes 9.
          REF_SIGN="intgr"; [[ "${mode}" == "sign_sepr" ]] && REF_SIGN="ss"
          REF_FIT="${MCDIR}dr_correction_fits_${LABEL}${WPS_SUF}${PTBIN_SUF}_step${step}_${method}$(sign_fsuf "${REF_SIGN}")${PMSUF}.root"
          NPT=$(fit_file_npt "${REF_FIT}" "${step}")
          EXP_PNG=$(( NPT + 1 ))
          # MAIN and RATIO canvases are counted SEPARATELY, and now live in two directories: the
          # main overlay canvases at the top of the view directory, the same-sign/opposite-sign
          # ratio canvases in its ratio/ subdir. Separately, because when one glob covered both
          # families it summed to 2*(NPT+1) against an EXP_PNG of NPT+1 -- so the check passed
          # with the ENTIRE ratio set missing, the one figure the sign split exists to produce.
          # sign_intgr has a single sign and no ratio canvas.
          # A MISSING directory means zero files, not a pipeline abort. `find` on a path that does
          # not exist returns 1, and under `set -o pipefail` that 1 survives the pipe to `wc` and
          # trips the ERR trap -- which is exactly how the first run of the two-view layout died
          # on the sign_intgr view, where ratio/ legitimately never exists. Count what is there.
          count_png() {   # $1 = directory, $2 = -name pattern, $3.. = extra find predicates
            local dir="$1"; shift
            [[ -d "${dir}" ]] || { echo 0; return 0; }
            find "${dir}" -maxdepth 1 "$@" 2>/dev/null | wc -l
          }
          NPNG=$(count_png "${MDIR}${mode}/${VIEWDIR}" -name "step${step}_dr_fit_${method}_*.png" \
                           ! -name "*_ratio.png")
          # MIRRORS THE C++ GATE (plot_dr_correction_fits.cxx: `rodir`): the ratio/ level is part
          # of the step-3 restructure only, so in step 4 the ratio canvases sit beside the main
          # ones. Checking for a ratio/ that the macro is not asked to write reported 6 x "0/9
          # ratio PNGs" for output that was in fact complete.
          RATIO_DIR="${MDIR}${mode}/${VIEWDIR}"
          [[ "${step}" == "3" ]] && RATIO_DIR="${MDIR}${mode}/${VIEWDIR}ratio"
          NPNG_RATIO=$(count_png "${RATIO_DIR}" \
                                 -name "step${step}_dr_fit_${method}_*_ratio.png")
          EXP_PNG_RATIO=0
          [[ "${mode}" == "sign_sepr" ]] && EXP_PNG_RATIO="${EXP_PNG}"
          if [[ "${NPT}" -lt 1 ]]; then
            ARTEFACT_FAILURES+=("cannot read the pair-pT binning from ${REF_FIT}")
          else
            if [[ "${NPNG}" -lt "${EXP_PNG}" ]]; then
              ARTEFACT_FAILURES+=("only ${NPNG}/${EXP_PNG} main PNGs in ${MDIR}${mode}/${VIEWDIR}")
              log "  !! only ${NPNG} main PNGs (expected ${EXP_PNG}: ${NPT} pair-pT bins + inclusive)"
            fi
            if [[ "${NPNG_RATIO}" -lt "${EXP_PNG_RATIO}" ]]; then
              ARTEFACT_FAILURES+=("only ${NPNG_RATIO}/${EXP_PNG_RATIO} same-sign/opposite-sign ratio PNGs in ${RATIO_DIR}")
              log "  !! only ${NPNG_RATIO} ratio PNGs (expected ${EXP_PNG_RATIO}: ${NPT} pair-pT bins + inclusive)"
            fi
          fi
        done   # view
        done   # mode

        # PERSISTENCE is a hard failure (the compiled-TF1 read-back trap); FLATNESS beyond Rp is
        # a property of the chosen function and is reported, not enforced -- `expo` is expected
        # to miss it, and that is exactly the information the method comparison needs.
        # One read-back report per SERIES; it stays in <method>/ next to the fit reports.
        for sgn in ${FITTED_SIGNS}; do
          RB="${MDIR}readback_check$(sign_rtag "${sgn}").txt"
          if [[ ! -f "${RB}" ]]; then
            ARTEFACT_FAILURES+=("no $(basename "${RB}") in ${MDIR}")
          else
            if ! grep -qE "^checked [0-9]+ functions, 0 FAILED persistence" "${RB}"; then
              ARTEFACT_FAILURES+=("READ-BACK PERSISTENCE FAILURE in ${RB}")
              log "  !! READ-BACK FAILURE -- persisted functions do not evaluate correctly"
            fi
            FLATLINE=$(grep "^flatness beyond Rp" "${RB}" || true)
            [[ -n "${FLATLINE}" ]] && log "  [$(sign_text "${sgn}")] ${FLATLINE}"
          fi
        done
      done
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
  echo "     See <plot base>/step<N>_dr_fit/plateau_corrected/plateau_guard_report.txt for the"
  echo "     offending cells (the no-plateau-correction mode has no fatal tier -- it applies no"
  echo "     plateau, so |plateau-1| cannot disqualify a cell there)."
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
