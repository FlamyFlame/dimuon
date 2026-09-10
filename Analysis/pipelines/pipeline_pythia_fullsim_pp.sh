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
# MC-TRIG-EFF (Stage 3 + Stage 10) is OFF by default: it is a heavier, separate deliverable and
# its final plots OVERWRITE the canonical pp_trigger_efficiency/mc_based/ dir (the pipeline backs
# that dir up first). The chain was given a pp_full sample knob (2026-07-20); the PbPb-union-weight
# physics question in docs/tracking/mc_trigger_efficiency.md is PbPb-only and does not affect pp
# (pp uses the 2mu4 PRODUCT weight). Enable with ENABLE_MC_TRIG_EFF=1.
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
    # ALRB / atlasLocalSetup.sh is NOT `set -e`/`set -u` safe: under either it aborts the whole
    # pipeline *silently* inside the setup (it even prints "Warning: -e is set ... may cause
    # issues"). Relax both across the setup, then restore the pipeline's strict flags.
    set +eu
    source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh" --quiet
    lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"
    set -eu
    command -v root >/dev/null || fail "ROOT not on PATH after lsetup"
}

# entries in a named tree (0 if missing) -- the validation primitive.
# ROOT frequently exits non-zero even on success; under `set -e`+`pipefail` that would kill a
# `var=$(tree_entries ...)` assignment SILENTLY (the `|| fail` sits on the next line). So the
# function swallows ROOT's exit status and always returns a number on stdout, exit 0.
tree_entries() {  # $1=file $2=tree
    local n
    n=$(root -l -b -q -e "
        TFile* f = TFile::Open(\"$1\");
        if (!f || f->IsZombie()) { printf(\"ENT 0\n\"); return; }
        TTree* t = (TTree*)f->Get(\"$2\");
        printf(\"ENT %lld\n\", t ? t->GetEntries() : 0LL);
    " 2>/dev/null | grep -oP 'ENT \K[0-9]+' | tail -1 || true)
    echo "${n:-0}"
}

# 1 if a named histogram exists in a file AND has entries, else 0. Same contract as
# tree_entries: swallows ROOT's exit status, always prints a number, always exits 0.
# This is the check that distinguishes "the event loop ran" from "the event loop threw and
# ROOT swallowed it, leaving a freshly-written near-empty file" -- a key COUNT cannot, because
# a partially-flushed file still has keys.
hist_filled() {  # $1=file $2=histogram
    local n
    n=$(root -l -b -q -e "
        TFile* f = TFile::Open(\"$1\");
        if (!f || f->IsZombie()) { printf(\"HF 0\n\"); return; }
        TH1* h = (TH1*)f->Get(\"$2\");
        printf(\"HF %d\n\", (h && h->GetEntries() > 0) ? 1 : 0);
    " 2>/dev/null | grep -oP 'HF \K[0-9]+' | tail -1 || true)
    echo "${n:-0}"
}

# Freshness stamps. A stage is NEVER validated by file-exists alone: Stages 1-3 and 5 back the
# OLD product up with `cp -a` and leave it at the nominal path, and a C++ throw inside an RDF
# event loop is caught by TRint so the job still exits 0. Existence + non-empty therefore passes
# on LAST production's file, and the whole chain would run on the pre-change 4.0 GeV / -1.30 gap
# population with no crash and no warning. The stamp makes "did this stage actually rewrite the
# file" answerable. (`pipeline_pp_crossx.sh` already does this; this pipeline did not.)
STAMP_DIR="$(mktemp -d)"
trap 'rm -rf "${STAMP_DIR}"' EXIT
stamp_now() { local s="${STAMP_DIR}/$1"; : > "$s"; echo "$s"; }
require_fresh() {  # $1=stamp $2..=files
    local stamp="$1"; shift
    local f
    for f in "$@"; do
        [[ -f "$f" ]] || fail "missing $f"
        [[ "$f" -nt "$stamp" ]] || fail "STALE: $(basename "$f") was NOT rewritten by the stage that just ran.
          It predates this stage, so it is the PREVIOUS production's file -- almost certainly a
          C++ exception thrown inside the ROOT/RDF event loop, which TRint catches so the job
          still exits 0. Do NOT let the chain continue: every downstream correction would be
          measured on the old selection. Check the stage log for 'Error'/'Runtime error'/'throw'."
    done
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
    # With ENABLE_MC_TRIG_EFF=1 Stage 3 ALSO rewrites the two mc_trig products, so they belong in
    # this loop -- without them the comment above was simply false for the configuration the
    # muon-pT 4.5 rerun uses.
    BACKUP_TARGETS=("${PAIR_FILE}" "${SINGLE_FILE}")
    if (( ENABLE_MC_TRIG_EFF )); then
        BACKUP_TARGETS+=("${SAMPLE_DIR}/muon_pairs_pythia_fullsim_pp24${CUT}_mc_trig${SFX}.root")
        BACKUP_TARGETS+=("${SAMPLE_DIR}/muon_pairs_pythia_fullsim_pp24${CUT}_single_muon_mc_trig${SFX}.root")
    fi
    for f in "${BACKUP_TARGETS[@]}"; do
        [[ -f "$f" ]] && { bak="${f%.root}.bak_$(date +%Y%m%d_%H%M%S).root"; cp -a "$f" "$bak"; log "  backed up $(basename "$f") -> $(basename "$bak")"; }
    done

    # --dry-run = a REAL smoke test: run every stage end to end, but cap the events per chain.
    # (Only the full-sample run scripts honour NEVENTS_MAX; the test sample is small already.)
    NEV_ENV=()
    if (( DRY_RUN )) && [[ "$SAMPLE" == "full" ]]; then
        NEV_ENV=(env "NEVENTS_MAX=${SMOKE_NEVENTS}")
        log "  SMOKE TEST: nevents_max=${SMOKE_NEVENTS} per (slice,beam) chain"
    fi

    NTP_STAMP="$(stamp_now ntp)"

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
# FRESHNESS FIRST -- see the require_fresh comment. Only meaningful if this run actually ran the
# NTP stages; under SKIP_NTP=1 the files are deliberately last run's and there is nothing to
# compare against.
if (( ! SKIP_NTP )); then
    require_fresh "${NTP_STAMP}" "${PAIR_FILE}" "${SINGLE_FILE}"
    log "  NTP outputs are freshly written  ✅"
fi
n_pair=$(tree_entries "${PAIR_FILE}" "muon_pair_tree_kin0_sign2") || n_pair=0
[[ "${n_pair:-0}" -gt 0 ]] || fail "muon_pair_tree_kin0_sign2 is EMPTY in ${PAIR_FILE}"
log "  pair tree kin0_sign2: ${n_pair} entries  ✅"
n_single=$(tree_entries "${SINGLE_FILE}" "muon_tree") || n_single=0
[[ "${n_single:-0}" -gt 0 ]] || fail "muon_tree is EMPTY in ${SINGLE_FILE} (Stage 2 silent throw)"
log "  single-muon tree: ${n_single} entries  ✅"

# If MC trig-eff is on, VALIDATE the mc_trig NTP too. A ROOT stage that THROWS still exits 0, so
# the NTP script's `|| fail` misses it (this is exactly how the store_mc_trigger multi-file bug
# silently skipped the mc_trig NTP on 2026-07-20). Check the artefact, not the exit code.
if (( ENABLE_MC_TRIG_EFF )); then
    MCTRIG_PAIR="${SAMPLE_DIR}/muon_pairs_pythia_fullsim_pp24${CUT}_mc_trig${SFX}.root"
    MCTRIG_SINGLE="${SAMPLE_DIR}/muon_pairs_pythia_fullsim_pp24${CUT}_single_muon_mc_trig${SFX}.root"
    [[ -f "$MCTRIG_PAIR" ]] || fail "mc_trig NTP missing: ${MCTRIG_PAIR} — Stage 3 likely threw but exited 0 (ROOT swallows C++ exceptions). Check the log for 'Runtime error'/'store_mc_trigger'."
    [[ -f "$MCTRIG_SINGLE" ]] || fail "mc_trig single-muon NTP missing: ${MCTRIG_SINGLE} — Stage 3 likely threw but exited 0."
    if (( ! SKIP_NTP )); then
        require_fresh "${NTP_STAMP}" "$MCTRIG_PAIR" "$MCTRIG_SINGLE"
        log "  mc_trig NTP outputs are freshly written  ✅"
    fi
    n_mct=$(tree_entries "$MCTRIG_PAIR" "muon_pair_tree_kin0_sign2") || n_mct=0
    [[ "${n_mct:-0}" -gt 0 ]] || fail "mc_trig NTP ${MCTRIG_PAIR} is EMPTY (Stage 3 silent throw)."
    log "  mc_trig pair tree kin0_sign2: ${n_mct} entries  ✅"
fi

# --- Stage 5: RDF histogram filling -------------------------------------------------------
if (( SKIP_RDF )); then
    log "[Stage 5] SKIPPED (SKIP_RDF=1)"
else
    [[ -f "${HIST_FILE}" ]] && { bak="${HIST_FILE%.root}.bak_$(date +%Y%m%d_%H%M%S).root"; cp -a "${HIST_FILE}" "$bak"; log "  backed up $(basename "${HIST_FILE}")"; }
    RDF_STAMP="$(stamp_now rdf)"
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
if (( ! SKIP_RDF )); then
    require_fresh "${RDF_STAMP}" "${HIST_FILE}"
    log "  histogram file is freshly written  ✅"
fi
# A key COUNT is not a validation: a throw part-way through the event loop still flushes some
# keys, and RECREATE on a file that then throws immediately leaves a near-empty one (the
# 851-byte, 0-key Pb+Pb corpse of 2026-09-06 is the extreme case, but the partial case passes a
# key count). Probe a histogram that only a completed SIGNAL-REGION fill can have filled.
REQUIRED_HIST="h_pair_pt_single_b_pass_tight"
[[ "$(hist_filled "${HIST_FILE}" "${REQUIRED_HIST}")" == "1" ]] \
    || fail "${HIST_FILE} has no filled ${REQUIRED_HIST} — the RDF event loop threw and ROOT swallowed it (exit code 0 is not evidence). Nothing downstream may consume this file."
nkeys=$(root -l -b -q -e "TFile* f=TFile::Open(\"${HIST_FILE}\"); printf(\"NK %d\n\", f&&!f->IsZombie()? f->GetListOfKeys()->GetSize():0);" 2>/dev/null | grep -oP 'NK \K[0-9]+' | tail -1 || true); nkeys=${nkeys:-0}
log "  ${nkeys} histogram keys, ${REQUIRED_HIST} filled  ✅"

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

    log "[Stage 8] single-muon reco efficiency + reco distributions"
    root -l -b <<ROOTEOF || fail "single-muon reco-eff plotting failed"
.L plot_single_muon_reco_effcy.cxx+
plot_single_muon_reco_effcy("pp", "default", ${TIGHT_CPP}, ${IS_TEST_CPP});
.q
ROOTEOF
    popd >/dev/null
    # reco-distr macro is a .C run INTERPRETED (not ACLiC); it lives at the Analysis root.
    #
    # OPT-IN SINCE 2026-09-09 -- this macro carries a RETIRED selection and is not maintained.
    # Its `sel_op` applies `pair_pt > 8` (the signal cut moved to ParamsSet::signal_pair_pt_min
    # = 9 GeV), the per-muon one-sided `q*eta < 2.2` (replaced by the fiducial gap windows on
    # 2026-08-17) and `dr > 0.05` (REMOVED from the analysis on 2026-06-22). Running it as part of
    # the nominal pipeline produced a figure describing a selection the analysis has not used for
    # months, beside figures that do use the current one. Set RUN_RECO_DISTR=1 to run it anyway.
    if [[ "${RUN_RECO_DISTR:-0}" == "1" ]]; then
        log "  ##### STALE SELECTION: plot_reco_distr_singleb_vs_op_pp24.C applies pair_pt > 8,"
        log "  ##### q*eta < 2.2 and dr > 0.05 -- none of which is the current signal region. #####"
        pushd "${ANALYSIS_DIR}" >/dev/null
        root -l -b -q "plot_reco_distr_singleb_vs_op_pp24.C(${IS_TEST_CPP}, ${TIGHT_CPP})" \
            || fail "reco-distr plotting failed"
        popd >/dev/null
    else
        log "  [skipped] plot_reco_distr_singleb_vs_op_pp24.C -- retired selection; RUN_RECO_DISTR=1 to run"
    fi

    log "[Stage 9] differential crossx per pT-hat slice + kn contributor table"
    pushd "${PY_PLOT_DIR}" >/dev/null
    root -l -b <<ROOTEOF || fail "crossx plotting failed"
.L plot_pythia_fullsim_kn_pt_crossx.cxx+
plot_pythia_fullsim_kn_pt_crossx(${IS_TEST_CPP}, ${TIGHT_CPP});
.q
ROOTEOF
    root -l -b <<ROOTEOF || fail "kn contributor table failed"
.L make_kn_contributor_table.cxx+
make_kn_contributor_table(${IS_TEST_CPP});
.q
ROOTEOF
    popd >/dev/null
fi

# --- Stage 10: MC-based trigger efficiency (OPTIONAL) -------------------------------------
# Chain: FillMCTrigEffHists(step1) -> FitMCSinglesEffcy -> FillMCTrigEffHists(do_step3) ->
# plot_mc_trig_eff. Intermediate hists/fits carry the sample label (pp24_full) so they never
# clobber the TEST-sample ones; the FINAL plots go to the canonical
# dimuon_data/plots/pp_trigger_efficiency/mc_based/ (the full sample SUPERSEDES the test-sample
# deliverable) -- so that directory is BACKED UP first.
# Overlay / r17663 trig-eff are NOT touched (this pipeline is pp only).
if (( ENABLE_MC_TRIG_EFF )); then
    log "[Stage 10] MC-based trigger efficiency chain"
    if [[ "$SAMPLE" == "full" ]]; then TRIG_SAMPLE="pp_full"; else TRIG_SAMPLE="pp"; fi

    PP_TRIG_PLOTS="${DATA_ROOT}/dimuon_data/plots/pp_trigger_efficiency/mc_based"
    if [[ "$SAMPLE" == "full" && -d "$PP_TRIG_PLOTS" ]]; then
        bak="${PP_TRIG_PLOTS}.bak_testsample_$(date +%Y%m%d_%H%M%S)"
        cp -a "$PP_TRIG_PLOTS" "$bak"
        log "  backed up TEST-sample pp trig-eff plots -> $(basename "$bak")"
    fi

    pushd "${RDF_DIR}" >/dev/null
    log "  [10a] FillMCTrigEffHists step1 (${TRIG_SAMPLE}, WP=$([[ $USE_TIGHT_WP == 1 ]] && echo tight || echo medium))"
    root -l -b <<ROOTEOF || fail "FillMCTrigEffHists step1 failed"
.L FillMCTrigEffHists.cxx+
FillMCTrigEffHists("${TRIG_SAMPLE}", false, ${TIGHT_CPP});
.q
ROOTEOF
    log "  [10b] FitMCSinglesEffcy (${TRIG_SAMPLE})"
    root -l -b <<ROOTEOF || fail "FitMCSinglesEffcy failed"
.L FitMCSinglesEffcy.cxx+
FitMCSinglesEffcy("${TRIG_SAMPLE}", ${TIGHT_CPP});
.q
ROOTEOF
    log "  [10c] FillMCTrigEffHists step3 (ε_ΔR inputs)"
    root -l -b <<ROOTEOF || fail "FillMCTrigEffHists step3 failed"
.L FillMCTrigEffHists.cxx+
FillMCTrigEffHists("${TRIG_SAMPLE}", true, ${TIGHT_CPP});
.q
ROOTEOF
    popd >/dev/null

    log "  [10d] plot_mc_trig_eff (${TRIG_SAMPLE})"
    pushd "${ANALYSIS_DIR}/plotting_codes/trig_effcy/mc_based" >/dev/null
    root -l -b <<ROOTEOF || fail "plot_mc_trig_eff failed"
.L plot_mc_trig_eff.cxx+
plot_mc_trig_eff("${TRIG_SAMPLE}", ${TIGHT_CPP});
.q
ROOTEOF
    popd >/dev/null
    log "[Stage 10] MC trig-eff done -> ${PP_TRIG_PLOTS}"
else
    log "[Stage 10] SKIPPED — MC trig-eff disabled (ENABLE_MC_TRIG_EFF=0; see the header note)"
fi

log "══════════ DONE — SAMPLE=${SAMPLE} ══════════"
log "  pair tree : ${PAIR_FILE}"
log "  hists     : ${HIST_FILE}"
log "  plots     : ${SAMPLE_DIR}/plots/"
