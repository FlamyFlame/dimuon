#!/usr/bin/env bash
# ============================================================================================
# RETIRED 2026-08-05 (advisor, relayed by the user): weighting each MC muon by SF = eps_data/
# eps_MC is NOT a reasonable procedure. The analysis uses data-driven single-muon efficiencies
# with MC-derived dR corrections; if the SF is ever to be tested it is applied as a direct
# multiplication eps_MC x SF, never as a per-muon weight.
# The R16 CONCLUSION stands (the dR corrections are insensitive to the single-muon
# normalization) -- only this implementation is retired. Its plot directories
# (step{1,3,4}_corrected_mc) have been deleted and are no longer produced.
# Kept executable only to reproduce the historical study on request.
# ============================================================================================
echo "run_mc_trigeff_corrected.sh is RETIRED (see header). Set I_MEAN_IT=1 to run anyway." >&2
[[ -z "${I_MEAN_IT:-}" ]] && exit 0
# =============================================================================
# CORRECTED-MC trigger-efficiency study driver
# (docs/tracking/mc_trigger_efficiency.md, round-7 Autonomy Contract item 5).
#
# Regenerates the COMPLETE corrected-MC product set from scratch. Question asked:
#   data and MC disagree on the single-muon trigger efficiency, so there are two ways to apply
#   the trigger correction --
#     (1) data-driven single-muon ε  +  MC-derived (ΔR, pair pT, pair η) corrections   [nominal]
#     (2) correct the MC to the data by weighting each MC muon by SF = ε_data/ε_MC(pT, q·η)
#   and we check (a) that the two give the SAME single-muon efficiency (Step 1) and (b) that
#   correcting the MC does NOT change the ΔR correction terms (Steps 3 and 4) -- the latter is
#   what makes the whole procedure valid, since the ΔR correction is measured in the MC world
#   and applied in the data world.
#
# Stages, per (sample, WP):
#   1. corrected Step-1 fill   -> mc_trig_eff_hists_<label><wp>_corrected.root
#   2. corrected turn-on fit   -> single_mu_effcy_pT_fit_mc_corrected<wp>.root
#   3. corrected Step-3 fill   -> mc_trig_eff_hists_<label><wp>_corrected_step3.root
#   4. corrected Step-4 fill   -> mc_trig_eff_hists_<label><wp>_corrected_step4.root
#   5. plots                   -> <plot base>/step{1,3,4}_corrected_mc/[medium/]
# Stage 3/4 read the CORRECTED fit of stage 2 (§3.3/§3.4 self-consistency), so the order is
# mandatory. NOTHING nominal is read for writing or overwritten: every artefact this script
# produces has "corrected" in its name.
#
# CLOSURE (`STAGES=closure`): reruns the same chain with SF forced to 1 (`sf_closure`), whose
# output must reproduce the NOMINAL fill bin-by-bin. It is the end-to-end test that the shared
# code path is unchanged; run it after ANY edit to FillMCTrigEffHists.cxx.
#
# Usage:  ./run_mc_trigeff_corrected.sh                     (everything, both samples, both WPs)
#         SAMPLES=overlay WPS=tight ./run_mc_trigeff_corrected.sh
#         STAGES="plots"   ./run_mc_trigeff_corrected.sh     (replot only)
#         STAGES="closure" ./run_mc_trigeff_corrected.sh
# =============================================================================
set -Euo pipefail

SAMPLES="${SAMPLES:-pp_full overlay}"
WPS="${WPS:-tight medium}"
STAGES="${STAGES:-fill plots}"          # fill | plots | closure

ANALYSIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
PLOT_DIR="${ANALYSIS_DIR}/plotting_codes/trig_effcy/mc_based"
LOG_DIR="${ANALYSIS_DIR}/pipelines/logs_corrected"
mkdir -p "$LOG_DIR"

ts(){ date "+%Y-%m-%d %H:%M:%S"; }
log(){ echo "[$(ts)] $*"; }
fail(){ echo "[$(ts)] FATAL: $*" >&2; exit 1; }
has_stage(){ [[ " $STAGES " == *" $1 "* ]]; }

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
set +eu
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh" --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"
set -eu

# sample -> MC output directory + hist-file label (mirrors FillMCTrigEffHists::GetSampleConfig)
outdir_of(){ case "$1" in
  pp_full) echo /usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample ;;
  overlay) echo /usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample ;;
  noovl)   echo /usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_no_overlay_test_sample ;;
  *) fail "unknown sample $1" ;; esac; }
label_of(){ case "$1" in
  pp_full) echo pp24_full ;;
  overlay) echo hijing_overlay_pbpb23 ;;
  noovl)   echo r17663_no_overlay ;;
  *) fail "unknown sample $1" ;; esac; }
# plot base (mirrors plot_mc_trig_eff_corrected.cxx MakeCfg)
plotbase_of(){ case "$1" in
  pp_full) echo /usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/mc_based ;;
  overlay) echo /usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pbpb_trigger_efficiency/mc_based ;;
  noovl)   echo /usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/r17663_no_overlay_trigger_efficiency/mc_based ;;
  *) fail "unknown sample $1" ;; esac; }

# A ROOT macro that throws still exits 0 -- NEVER trust the exit code, always check the artefact.
val_file(){ [[ -s "$1" ]] || fail "missing/empty artefact: $1"; }

# ---- 0. pre-compile the ACLiC macros ONCE (serial: a race on the shared _cxx.so is the one
#         real hazard of running the per-sample chains in parallel) -----------------------------
log "════ pre-compiling ACLiC macros ════"
# NB: `root -q 'X.cxx+'` compiles AND RUNS X() with default arguments -- use `.L X.cxx+`.
( cd "$RDF_DIR"  && root -l -b -q -e '.L FillMCTrigEffHists.cxx+' ) >"$LOG_DIR/compile_fill.log" 2>&1 \
    || fail "FillMCTrigEffHists compile failed (see $LOG_DIR/compile_fill.log)"
( cd "$RDF_DIR"  && root -l -b -q -e '.L FitMCSinglesEffcy.cxx+'  ) >"$LOG_DIR/compile_fit.log"  2>&1 \
    || fail "FitMCSinglesEffcy compile failed (see $LOG_DIR/compile_fit.log)"
( cd "$PLOT_DIR" && root -l -b -q -e '.L plot_mc_trig_eff_corrected.cxx+' ) \
    >"$LOG_DIR/compile_plot.log" 2>&1 || fail "plot_mc_trig_eff_corrected compile failed"
log "  compiled OK"

# ---- 1. corrected fill / fit chain, one background job per (sample, WP) ------------------------
chain(){   # $1=sample  $2=wp  $3=closure(true|false)
    local s="$1" wp="$2" clo="$3" cpp suf out lbl lg tag
    cpp=$([[ $wp == tight ]] && echo true || echo false)
    suf=$([[ $wp == tight ]] && echo ""   || echo "_medium_wp")
    tag=$([[ $clo == true ]] && echo "_corrected_sfclosure" || echo "_corrected")
    out="$(outdir_of "$s")"; lbl="$(label_of "$s")"
    lg="$LOG_DIR/${s}_${wp}${tag}.log"
    {
        echo "=== $s / $wp / closure=$clo : corrected step1 ==="
        ( cd "$RDF_DIR" && root -l -b -q \
            "FillMCTrigEffHists.cxx+(\"$s\", false, $cpp, false, false, true, $clo)" )
        val_file "$out/mc_trig_eff_hists_${lbl}${suf}${tag}.root"
        echo "=== $s / $wp / closure=$clo : corrected fit ==="
        # The fit output goes to its OWN file first, then is echoed into this log and grepped.
        # (`| tee /dev/stderr` would re-open the log via /proc/self/fd/2 at offset 0 and overwrite
        #  everything written so far -- it silently truncated these logs once already.)
        local fitlog="$LOG_DIR/fit_${s}_${wp}${tag}.log"
        ( cd "$RDF_DIR" && root -l -b -q "FitMCSinglesEffcy.cxx+(\"$s\", $cpp, true, $clo)" ) \
            >"$fitlog" 2>&1
        cat "$fitlog"
        # An EMPTY fit graph leaves the TF1 at its initial parameters while still returning
        # status 0, so a silent bad turn-on must be caught here, not downstream.
        grep -q "fits, 0 failed" "$fitlog" \
            || fail "corrected turn-on fit reported FAILED fits for $s/$wp (see $fitlog)"
        val_file "$out/single_mu_effcy_pT_fit_mc${tag}${suf}.root"
        echo "=== $s / $wp / closure=$clo : corrected step3 ==="
        ( cd "$RDF_DIR" && root -l -b -q \
            "FillMCTrigEffHists.cxx+(\"$s\", true, $cpp, false, false, true, $clo)" )
        val_file "$out/mc_trig_eff_hists_${lbl}${suf}${tag}_step3.root"
        echo "=== $s / $wp / closure=$clo : corrected step4 ==="
        ( cd "$RDF_DIR" && root -l -b -q \
            "FillMCTrigEffHists.cxx+(\"$s\", false, $cpp, true, false, true, $clo)" )
        val_file "$out/mc_trig_eff_hists_${lbl}${suf}${tag}_step4.root"
        echo "=== $s / $wp / closure=$clo : CHAIN OK ==="
    } >"$lg" 2>&1
}

run_chains(){   # $1 = closure flag
    local clo="$1" rc=0
    local pids=() names=()
    for s in $SAMPLES; do for wp in $WPS; do
        log "launching corrected chain (closure=$clo): $s / $wp"
        chain "$s" "$wp" "$clo" & pids+=($!); names+=("$s/$wp")
    done; done
    for i in "${!pids[@]}"; do
        if wait "${pids[$i]}"; then log "  chain OK : ${names[$i]}"
        else log "  chain FAILED : ${names[$i]}"; rc=1; fi
    done
    [[ $rc -eq 0 ]] || fail "at least one corrected chain failed -- see $LOG_DIR/*.log"
}

if has_stage fill;    then log "════ corrected fill/fit chains ════"; run_chains false; fi
if has_stage closure; then log "════ SF≡1 CLOSURE chains ════";       run_chains true;  fi

# ---- 2. plots -------------------------------------------------------------------------------
if has_stage plots; then
  log "════ corrected plots ════"
  for wp in $WPS; do
    cpp=$([[ $wp == tight ]] && echo true || echo false)
    wpdir=$([[ $wp == tight ]] && echo "" || echo "medium/")
    for s in $SAMPLES; do
        log "plotting $s / $wp"
        ( cd "$PLOT_DIR" && root -l -b -q "plot_mc_trig_eff_corrected.cxx+(\"$s\", $cpp)" ) \
            >"$LOG_DIR/plot_${s}_${wp}.log" 2>&1 || fail "plot $s/$wp failed"
        grep -q "plot_mc_trig_eff_corrected(.*) done\." "$LOG_DIR/plot_${s}_${wp}.log" \
            || fail "plot $s/$wp did not reach the end (see $LOG_DIR/plot_${s}_${wp}.log)"
        # artefact validation (a throwing macro still exits 0)
        pb="$(plotbase_of "$s")"
        for f in "step1_corrected_mc/${wpdir}step1_corrected_eff_pt.png" \
                 "step1_corrected_mc/${wpdir}step1_corrected_vs_data.txt" \
                 "step3_corrected_mc/${wpdir}step3_corr_vs_orig_zoom_pairpt1.png" \
                 "step3_corrected_mc/${wpdir}step3_corr_vs_orig_pulls.txt" \
                 "step4_corrected_mc/${wpdir}step4_corr_vs_orig_zoom_pairpt1.png" \
                 "step4_corrected_mc/${wpdir}step4_corr_vs_orig_pulls.txt" ; do
            val_file "$pb/$f"
        done
    done
  done
fi

log "════ corrected-MC study DONE ($STAGES) ════"
