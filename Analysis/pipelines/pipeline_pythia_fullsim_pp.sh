#!/usr/bin/env bash
set -Eeuo pipefail

# ==========================================================================================
# End-to-end Pythia fullsim **pp24-conditions** pipeline (there was none before this).
#
#   Stage 1  NTuple processing  -- nominal pair tree
#   Stage 2  NTuple processing  -- single-muon tree
#   Stage 3  NTuple processing  -- MC-trigger variants (pair + single-muon)   [OPTIONAL]
#   Stage 4  Validate NTP outputs
#   Stage 5  RDF histogram filling (reco-eff + detector-response hists)
#   Stage 6  Validate histogram output
#   Stage 7  Reco efficiency + detector response plots
#   Stage 8  Single-muon reco efficiency plots
#   Stage 9  Differential crossx per pT-hat slice (the "statistics" plot) + kn contributor table
#   Stage 10 MC-based trigger efficiency chain                                 [OPTIONAL, default OFF]
#
# SAMPLE selection (the single switch that also fixes the ISOSPIN treatment -- see
# MuonObjectsParamsAndHelpers/FullSimSampleType.h and Analysis/docs/ami_weights.md):
#   full : the "_pdf" production, DSIDs 803015-803020, read from the LOCALGROUPDISK symlink farm.
#          pp beam ONLY, isospin weight 1  =>  an HONEST pp cross-section.
#          Products end in "_full"; everything is written under pythia_fullsim_full_sample/.
#   test : the small 4-beam TEST sample (produced with 4 isospin beams BY MISTAKE, so its absolute
#          sigma carries the Pb 4:6:6:9 isospin AVERAGE and is NOT a physical pp cross-section).
#
# Usage:
#   ./pipeline_pythia_fullsim_pp.sh full
#   ./pipeline_pythia_fullsim_pp.sh test
#   ./pipeline_pythia_fullsim_pp.sh full --dry-run          # tiny nevents_max smoke test
#   ENABLE_MC_TRIG_EFF=1 ./pipeline_pythia_fullsim_pp.sh full
#
# Optional env vars:
#   ENABLE_MC_TRIG_EFF=0   # Stage 3 + Stage 10. DEFAULT OFF -- see the note below.
#   USE_TIGHT_WP=1         # NOMINAL muon WP = Tight. 0 => Medium (distinct output dirs).
#   SMOKE_NEVENTS=2000     # events per (slice,beam) chain when --dry-run
#   SKIP_NTP=0 SKIP_RDF=0 SKIP_PLOTS=0
#
# ⚠ WHY MC-TRIG-EFF IS OFF BY DEFAULT (2026-07-14): the MC-trigger-efficiency chain
# (FillMCTrigEffHists / FitMCSinglesEffcy / plot_mc_trig_eff) is owned by a concurrently-running
# effort tracked in docs/tracking/mc_trigger_efficiency.md, which has an OPEN physics decision (the
# single-muon efficiency is DeltaR-dependent => the PbPb mu4 UNION weight is at risk) and whose
# sample directory is still hard-coded to the TEST sample. Turn this on only once that work has
# landed and the chain has been given a sample knob.
# ==========================================================================================

SAMPLE="${1:-}"
DRY_RUN=0
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
    esac
done

case "$SAMPLE" in
    full|test) ;;
    *) echo "Usage: $0 [full|test] [--dry-run]"; exit 1 ;;
esac

ENABLE_MC_TRIG_EFF="${ENABLE_MC_TRIG_EFF:-0}"
USE_TIGHT_WP="${USE_TIGHT_WP:-1}"
SMOKE_NEVENTS="${SMOKE_NEVENTS:-2000}"
SKIP_NTP="${SKIP_NTP:-0}"
SKIP_RDF="${SKIP_RDF:-0}"
SKIP_PLOTS="${SKIP_PLOTS:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
NTP_DIR="${ANALYSIS_DIR}/NTupleProcessingCode"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
RECO_DIR="${ANALYSIS_DIR}/plotting_codes/reco_effcy"
PY_PLOT_DIR="${ANALYSIS_DIR}/plotting_codes/pythia_plotting_codes"

DATA_ROOT="/usatlas/u/yuhanguo/usatlasdata"
if [[ "$SAMPLE" == "full" ]]; then
    SAMPLE_DIR="${DATA_ROOT}/pythia_fullsim_full_sample"
    SFX="_full"
    IS_TEST_CPP="false"
else
    SAMPLE_DIR="${DATA_ROOT}/pythia_fullsim_test_sample"
    SFX=""
    IS_TEST_CPP="true"
fi
CUT="_no_data_resonance_cuts"
PAIR_FILE="${SAMPLE_DIR}/muon_pairs_pythia_fullsim_pp24${CUT}${SFX}.root"
SINGLE_FILE="${SAMPLE_DIR}/muon_pairs_pythia_fullsim_pp24${CUT}_single_muon${SFX}.root"
HIST_FILE="${SAMPLE_DIR}/histograms_pythia_fullsim_pp24${CUT}${SFX}.root"

ts()   { date "+%Y-%m-%d %H:%M:%S"; }
log()  { echo "[$(ts)] $*"; }
fail() { echo "[$(ts)] FATAL: $*" >&2; exit 1; }

setup_root() {
    export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
    set +u
    source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh" --quiet
    lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"
    set -u
}

# entries in a named tree (0 if missing) -- the validation primitive
tree_entries() {  # $1=file $2=tree
    root -l -b -q -e "
        TFile* f = TFile::Open(\"$1\");
        if (!f || f->IsZombie()) { printf(\"ENT 0\n\"); return; }
        TTree* t = (TTree*)f->Get(\"$2\");
        printf(\"ENT %lld\n\", t ? t->GetEntries() : 0LL);
    " 2>/dev/null | grep -oP 'ENT \K[0-9]+' | tail -1
}

log "══════════ pp24 Pythia fullsim pipeline — SAMPLE=${SAMPLE} ══════════"
log "  sample dir      : ${SAMPLE_DIR}"
log "  WP              : $([[ $USE_TIGHT_WP == 1 ]] && echo Tight || echo Medium)"
log "  MC trig-eff     : $([[ $ENABLE_MC_TRIG_EFF == 1 ]] && echo ON || echo OFF)"
log "  dry-run (smoke) : ${DRY_RUN}$([[ $DRY_RUN == 1 ]] && echo "  (nevents_max=${SMOKE_NEVENTS})")"

# --- Stage 0: preflight --------------------------------------------------------------------
if [[ "$SAMPLE" == "full" ]]; then
    # ALL SIX pT-hat slices must be present. A missing slice is a BIAS, not a statistics loss:
    # it drops an entire sigma-weighted term from the combination. (The NTP code also throws --
    # allow_missing_slices=false -- but fail here, earlier and with a clearer message.)
    missing=()
    for s in pTH8_14 pTH14_24 pTH24_40 pTH40_70 pTH70_125 pTH125_300; do
        ls "${SAMPLE_DIR}/Pythia_5p36TeV_pp_hQCD_DiMu_${s}.FullSimPP24.NTUP.part"*.root >/dev/null 2>&1 \
            || missing+=("$s")
    done
    if (( ${#missing[@]} )); then
        fail "LOCALGROUPDISK farm INCOMPLETE — ${#missing[@]}/6 pT-hat slices missing: ${missing[*]}
          A missing slice BIASES the cross-section-weighted combination (it is not merely a
          statistics loss). Wait for the grid tasks + SkimCode/scripts/fullsim_pp24_full_to_lgd.sh."
    fi
    n_ami=$(ls "${SAMPLE_DIR}"/ami_info/ami_info_*.txt 2>/dev/null | wc -l)
    [[ "$n_ami" -eq 6 ]] || fail "expected 6 AMI files in ${SAMPLE_DIR}/ami_info, found ${n_ami} (see docs/ami_weights.md — a WRONG or MISSING AMI weight is a silent, non-cancelling error)"
    log "[Stage 0] preflight OK: all 6 farm slices present, 6 AMI files present"
fi

setup_root

# --- Stages 1-3: NTuple processing --------------------------------------------------------
if (( SKIP_NTP )); then
    log "[Stage 1-3] SKIPPED (SKIP_NTP=1)"
else
    if [[ "$SAMPLE" == "full" ]]; then
        NTP_NOMINAL="run_pythia_fullsim_full_sample.sh"
        NTP_SINGLE="run_pythia_fullsim_single_muon_full_sample.sh"
        NTP_TRIG="run_pythia_fullsim_mc_trig_full_sample.sh"
        NTP_TRIG_SINGLE="run_pythia_fullsim_single_muon_mc_trig_full_sample.sh"
    else
        NTP_NOMINAL="run_pythia_fullsim.sh"
        NTP_SINGLE="run_pythia_fullsim_single_muon.sh"
        NTP_TRIG="run_pythia_fullsim_mc_trig.sh"
        NTP_TRIG_SINGLE="run_pythia_fullsim_single_muon_mc_trig.sh"
    fi

    # Back up anything we are about to overwrite. Never clobber a previous result silently.
    for f in "${PAIR_FILE}" "${SINGLE_FILE}"; do
        [[ -f "$f" ]] && { bak="${f%.root}.bak_$(date +%Y%m%d_%H%M%S).root"; cp -a "$f" "$bak"; log "  backed up $(basename "$f") -> $(basename "$bak")"; }
    done

    # --dry-run = a REAL smoke test: run every stage end to end, but cap the events per chain.
    # (Only the full-sample run scripts honour NEVENTS_MAX; the test sample is small already.)
    NEV_ENV=()
    if (( DRY_RUN )) && [[ "$SAMPLE" == "full" ]]; then
        NEV_ENV=(env "NEVENTS_MAX=${SMOKE_NEVENTS}")
        log "  SMOKE TEST: nevents_max=${SMOKE_NEVENTS} per (slice,beam) chain"
    fi

    log "[Stage 1] NTP nominal   (${NTP_NOMINAL})"
    ( cd "${NTP_DIR}" && "${NEV_ENV[@]}" bash "${NTP_NOMINAL}" ) || fail "NTP nominal failed"

    log "[Stage 2] NTP single-mu (${NTP_SINGLE})"
    ( cd "${NTP_DIR}" && "${NEV_ENV[@]}" bash "${NTP_SINGLE}" ) || fail "NTP single-muon failed"

    if (( ENABLE_MC_TRIG_EFF )); then
        log "[Stage 3] NTP mc-trig  (${NTP_TRIG}, ${NTP_TRIG_SINGLE})"
        ( cd "${NTP_DIR}" && "${NEV_ENV[@]}" bash "${NTP_TRIG}" )        || fail "NTP mc_trig failed"
        ( cd "${NTP_DIR}" && "${NEV_ENV[@]}" bash "${NTP_TRIG_SINGLE}" ) || fail "NTP mc_trig single-muon failed"
    else
        log "[Stage 3] SKIPPED — MC trig-eff disabled (ENABLE_MC_TRIG_EFF=0)"
    fi
fi

# --- Stage 4: validate NTP outputs --------------------------------------------------------
log "[Stage 4] validating NTP output"
[[ -f "${PAIR_FILE}" ]]   || fail "missing ${PAIR_FILE}"
[[ -f "${SINGLE_FILE}" ]] || fail "missing ${SINGLE_FILE}"
n_pair=$(tree_entries "${PAIR_FILE}" "muon_pair_tree_kin0_sign2")
[[ "${n_pair:-0}" -gt 0 ]] || fail "muon_pair_tree_kin0_sign2 is EMPTY in ${PAIR_FILE}"
log "  pair tree kin0_sign2: ${n_pair} entries  ✅"

# --- Stage 5: RDF histogram filling -------------------------------------------------------
if (( SKIP_RDF )); then
    log "[Stage 5] SKIPPED (SKIP_RDF=1)"
else
    [[ -f "${HIST_FILE}" ]] && { bak="${HIST_FILE%.root}.bak_$(date +%Y%m%d_%H%M%S).root"; cp -a "${HIST_FILE}" "$bak"; log "  backed up $(basename "${HIST_FILE}")"; }
    log "[Stage 5] RDF histogram filling"
    pushd "${RDF_DIR}" >/dev/null
    root -l -b <<ROOTEOF || fail "ROOT exited non-zero during RDF hist filling"
.L RDFBasedHistFillingPythiaFullsim.cxx+
{
    RDFBasedHistFillingPythiaFullsim fs;
    fs.is_test_sample = ${IS_TEST_CPP};   // drives input dir AND the "_full" file suffix
    fs.Run();
}
.q
ROOTEOF
    popd >/dev/null
fi

# --- Stage 6: validate histogram output ---------------------------------------------------
log "[Stage 6] validating histogram output"
[[ -f "${HIST_FILE}" ]] || fail "missing ${HIST_FILE}"
nkeys=$(root -l -b -q -e "TFile* f=TFile::Open(\"${HIST_FILE}\"); printf(\"NK %d\n\", f&&!f->IsZombie()? f->GetListOfKeys()->GetSize():0);" 2>/dev/null | grep -oP 'NK \K[0-9]+' | tail -1)
[[ "${nkeys:-0}" -gt 0 ]] || fail "${HIST_FILE} has no keys"
log "  ${nkeys} histogram keys  ✅"

# --- Stages 7-9: plots --------------------------------------------------------------------
if (( SKIP_PLOTS )); then
    log "[Stage 7-9] SKIPPED (SKIP_PLOTS=1)"
else
    TIGHT_CPP=$([[ $USE_TIGHT_WP == 1 ]] && echo true || echo false)

    log "[Stage 7] reco efficiency + detector response"
    pushd "${RECO_DIR}" >/dev/null
    root -l -b <<ROOTEOF || fail "reco-eff / det-response plotting failed"
.L PythiaFullsimRecoEffPlotter.cxx+
{
    PythiaFullsimRecoEffPlotter pl(${TIGHT_CPP});
    pl.is_test_sample = ${IS_TEST_CPP};
    pl.Run();
}
.q
ROOTEOF

    log "[Stage 8] single-muon reco efficiency"
    root -l -b <<ROOTEOF || fail "single-muon reco-eff plotting failed"
.L plot_single_muon_reco_effcy.cxx+
plot_single_muon_reco_effcy("pp", "default", ${TIGHT_CPP}, ${IS_TEST_CPP});
.q
ROOTEOF
    popd >/dev/null

    log "[Stage 9] differential crossx per pT-hat slice + kn contributor table"
    pushd "${PY_PLOT_DIR}" >/dev/null
    root -l -b <<ROOTEOF || fail "crossx plotting failed"
.L plot_pythia_fullsim_kn_pt_crossx.cxx+
g_is_test_sample = ${IS_TEST_CPP};
g_use_tight_wp   = ${TIGHT_CPP};
plot_pythia_fullsim_kn_pt_crossx();
.q
ROOTEOF
    popd >/dev/null
fi

# --- Stage 10: MC-based trigger efficiency (OPTIONAL) -------------------------------------
if (( ENABLE_MC_TRIG_EFF )); then
    log "[Stage 10] MC-based trigger efficiency chain"
    fail "Stage 10 is not wired yet: FillMCTrigEffHists / FitMCSinglesEffcy / plot_mc_trig_eff still
          hard-code the TEST-sample directory and have no sample knob. That chain is owned by the
          concurrent effort in docs/tracking/mc_trigger_efficiency.md (open physics decision). Add
          the sample knob once that work lands, then enable this stage."
else
    log "[Stage 10] SKIPPED — MC trig-eff disabled (see the header note)"
fi

log "══════════ DONE — SAMPLE=${SAMPLE} ══════════"
log "  pair tree : ${PAIR_FILE}"
log "  hists     : ${HIST_FILE}"
log "  plots     : ${SAMPLE_DIR}/plots/"
