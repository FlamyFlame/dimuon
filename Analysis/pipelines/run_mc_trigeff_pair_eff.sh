#!/usr/bin/env bash
# =================================================================================================
# SINGLE-VALUE PAIR 2mu4 EFFICIENCY -- measurement, closure, comparison figure
# (docs/tracking/mc_trigeff_single_value_pair_eff.md)
#
# Measures ONE 2mu4 efficiency per (pair pT coarse bin x |eta^pair| group x pair sign) cell inside
# a dimuon mass window, as the ALTERNATIVE to the factorized
#     eps(pT1,q*eta1) eps(pT2,q*eta2) eps_dR(dR ; cell)
# weight in the three highest pair-pT cells, where the dR-shape fit runs out of pairs. Two numbers
# per cell are delivered: eps^pair (PURE, replaces the whole weight) and K (CALIBRATED, multiplies
# the two single-muon efficiencies -- the single-number analogue of eps_dR).
#
# STAGES
#   Stage 0  compile the three macros ONCE (a race on a shared _cxx.so is the real hazard)
#   Stage 1  FillMCTrigEffPairEff   -> pair_trig_eff_<label><wp>.root      (the deliverable)
#   Stage 2  FillMCTrigEffClosure   -> mc_trig_eff_closure_...<mode>.root, for the TWO dR
#            approaches this doc compares against, now carrying the single-value numerators too
#   Stage 3  plot_mc_trig_eff_closure_highpt_compare
#            -> <plot base>/closure/single_value_highpt_comparison/<comparison>/*.png
#               two comparisons x two applied forms = 4 PNGs per WP:
#                 mass_window_compr/  the two mass windows on the canonical cells
#                 pt_merge_compr/     the signal window, with and without the pair-pT merge
#   Stage 4  write_pair_trig_eff_tables -> <plot base>/single_value_pair_eff_tables/*.csv
#
# UPSTREAM IS NOT REBUILT HERE. The Step-3 dR fits and the MC single-muon turn-ons come from the
# trigger-efficiency chain (pipeline_pythia_fullsim_pp.sh Stage 10 + run_dr_correction_fits.sh) and
# are used as they are; Stage 2 THROWS if they describe a different pair-eta axis from today's.
#
# WHY ARTEFACT VALIDATION AND NOT EXIT CODES: a ROOT macro that throws inside an RDF event loop
# still exits 0, leaving a fresh but near-empty output. Every stage is checked by what it produced.
#
# Usage:
#   ./run_mc_trigeff_pair_eff.sh                 # Tight (nominal)
#   WPS="tight medium" ./run_mc_trigeff_pair_eff.sh
#   SKIP_FILL=1 ./run_mc_trigeff_pair_eff.sh     # reuse the ROOT files, remake the figure only
#
# Env vars:
#   WPS="tight"        working points (registry: Analysis/docs/muon_wp_registry.md). Default Tight.
#   MODES="..."        the dR approaches the figure compares against. Both are required by the
#                      comparison macro; changing this is for debugging only.
#   SKIP_FILL=1        skip Stages 1-2.
# =================================================================================================
set -Eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
PLOT_DIR="${ANALYSIS_DIR}/plotting_codes/trig_effcy/mc_based"
LOG_DIR="${SCRIPT_DIR}/logs_pair_eff"
MC_DIR="/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample"
PLOT_BASE="/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency"
mkdir -p "${LOG_DIR}"

SAMPLE="pp_full"
LABEL="pp24_full"
WPS="${WPS:-tight}"
MODES="${MODES:-nocorr_ptmerge nocorr_etamerge_ptmerge}"
SKIP_FILL="${SKIP_FILL:-0}"

log()  { echo "[$(date '+%H:%M:%S')] $*"; }
fail() { echo "[$(date '+%H:%M:%S')] FAIL: $*" >&2; exit 1; }

# A ROOT file is only accepted when it is NEWER than the run started AND carries a required key --
# both, because a throw leaves a fresh file behind and a reused file carries the key.
check_artifact() {  # <file> <key> <start-epoch> <what>
    local f="$1" key="$2" t0="$3" what="$4"
    [[ -f "$f" ]] || fail "${what}: ${f} was not produced"
    local mt; mt=$(stat -c %Y "$f")
    (( mt >= t0 )) || fail "${what}: ${f} is older than this run (mtime ${mt} < ${t0}) -- the macro threw"
    root -l -b -q -e "TFile* f = TFile::Open(\"${f}\"); if (!f || f->IsZombie() || !f->Get(\"${key}\")) { printf(\"MISSINGKEY\\n\"); } else printf(\"KEYOK\\n\");" 2>/dev/null \
        | grep -q KEYOK || fail "${what}: ${f} has no '${key}' -- the macro threw inside the event loop"
}

wp_suffix() { [[ "$1" == "tight" ]] && echo "" || echo "_medium_wp"; }
wp_cpp()    { [[ "$1" == "tight" ]] && echo "true" || echo "false"; }

mode_tag() {  # the DrCorrPlateauModeTag convention: "_" + token
    echo "_$1"
}

log "=== single-value pair 2mu4 efficiency ==="
log "  sample ${SAMPLE}, WPs '${WPS}', dR approaches '${MODES}'"

# --- Stage 0: compile ---------------------------------------------------------------------------
log "[Stage 0] compiling"
( cd "${RDF_DIR}"  && root -l -b -q -e '.L FillMCTrigEffPairEff.cxx+' ) \
    > "${LOG_DIR}/compile_paireff.log" 2>&1 || fail "FillMCTrigEffPairEff did not compile"
( cd "${RDF_DIR}"  && root -l -b -q -e '.L FillMCTrigEffClosure.cxx+' ) \
    > "${LOG_DIR}/compile_closure.log" 2>&1 || fail "FillMCTrigEffClosure did not compile"
( cd "${PLOT_DIR}" && root -l -b -q -e '.L plot_mc_trig_eff_closure_highpt_compare.cxx+' ) \
    > "${LOG_DIR}/compile_plot.log" 2>&1 || fail "plot_mc_trig_eff_closure_highpt_compare did not compile"
( cd "${PLOT_DIR}" && root -l -b -q -e '.L write_pair_trig_eff_tables.cxx+' ) \
    > "${LOG_DIR}/compile_tables.log" 2>&1 || fail "write_pair_trig_eff_tables did not compile"

for WP in ${WPS}; do
    SUF=$(wp_suffix "$WP"); CPP=$(wp_cpp "$WP")

    if (( ! SKIP_FILL )); then
        # --- Stage 1: the measurement ------------------------------------------------------------
        T0=$(date +%s)
        log "[Stage 1] FillMCTrigEffPairEff (${SAMPLE}, ${WP})"
        ( cd "${RDF_DIR}" && root -l -b <<ROOTEOF
.L FillMCTrigEffPairEff.cxx+
FillMCTrigEffPairEff("${SAMPLE}", ${CPP});
.q
ROOTEOF
        ) > "${LOG_DIR}/paireff_${WP}.log" 2>&1 || fail "FillMCTrigEffPairEff (${WP}) crashed"
        check_artifact "${MC_DIR}/pair_trig_eff_${LABEL}${SUF}.root" "h_paireff_eps_os_sig" \
                       "$T0" "Stage 1 (${WP})"
        log "  -> ${MC_DIR}/pair_trig_eff_${LABEL}${SUF}.root"

        # --- Stage 2: the closure, one run per dR approach ---------------------------------------
        for MODE in ${MODES}; do
            T0=$(date +%s)
            log "[Stage 2] FillMCTrigEffClosure (${SAMPLE}, ${WP}, mode ${MODE})"
            ( cd "${RDF_DIR}" && root -l -b <<ROOTEOF
.L FillMCTrigEffClosure.cxx+
FillMCTrigEffClosure("${SAMPLE}", ${CPP}, "${MODE}");
.q
ROOTEOF
            ) > "${LOG_DIR}/closure_${WP}_${MODE}.log" 2>&1 \
                || fail "FillMCTrigEffClosure (${WP}, ${MODE}) crashed"
            check_artifact "${MC_DIR}/mc_trig_eff_closure_${LABEL}${SUF}$(mode_tag "$MODE").root" \
                           "h_closure_signal_num_paireff_sig" "$T0" "Stage 2 (${WP}, ${MODE})"
        done
    else
        log "[Stages 1-2] SKIPPED (SKIP_FILL=1)"
    fi

    # --- Stage 3: the comparison figure ----------------------------------------------------------
    T0=$(date +%s)
    log "[Stage 3] plot_mc_trig_eff_closure_highpt_compare (${WP})"
    ( cd "${PLOT_DIR}" && root -l -b <<ROOTEOF
.L plot_mc_trig_eff_closure_highpt_compare.cxx+
plot_mc_trig_eff_closure_highpt_compare("${SAMPLE}", ${CPP});
.q
ROOTEOF
    ) > "${LOG_DIR}/plot_${WP}.log" 2>&1 || fail "the comparison figure (${WP}) crashed"
    BASE="${PLOT_BASE}/mc_based$( [[ "$WP" == medium ]] && echo _medium )"
    OUT="${BASE}/closure/single_value_highpt_comparison"
    for C in mass_window_compr pt_merge_compr; do
        for F in "closure_highpt_single_value_${C}_pure.png" \
                 "closure_highpt_single_value_${C}_calibrated.png"; do
            [[ -f "${OUT}/${C}/${F}" ]] || fail "Stage 3 (${WP}): ${OUT}/${C}/${F} was not written"
            mt=$(stat -c %Y "${OUT}/${C}/${F}")
            (( mt >= T0 )) || fail "Stage 3 (${WP}): ${C}/${F} is stale"
        done
    done
    log "  -> ${OUT}/{mass_window_compr,pt_merge_compr}/"

    # --- Stage 4: the CSV tables ------------------------------------------------------------------
    T0=$(date +%s)
    log "[Stage 4] write_pair_trig_eff_tables (${WP})"
    ( cd "${PLOT_DIR}" && root -l -b <<ROOTEOF
.L write_pair_trig_eff_tables.cxx+
write_pair_trig_eff_tables("${SAMPLE}", ${CPP});
.q
ROOTEOF
    ) > "${LOG_DIR}/tables_${WP}.log" 2>&1 || fail "the CSV tables (${WP}) crashed"
    TBL="${BASE}/single_value_pair_eff_tables"
    for F in single_value_pair_eff_opposite_sign${SUF}.csv \
             single_value_pair_eff_same_sign${SUF}.csv \
             single_value_pair_eff_opposite_sign_ptmerge${SUF}.csv \
             single_value_pair_eff_same_sign_ptmerge${SUF}.csv \
             single_value_pair_stats_same_sign${SUF}.csv \
             single_value_pair_stats_same_sign_ptmerge${SUF}.csv; do
        [[ -f "${TBL}/${F}" ]] || fail "Stage 4 (${WP}): ${TBL}/${F} was not written"
        mt=$(stat -c %Y "${TBL}/${F}"); (( mt >= T0 )) || fail "Stage 4 (${WP}): ${F} is stale"
    done
    log "  -> ${TBL}/"
done

log "=== done. logs in ${LOG_DIR} ==="
