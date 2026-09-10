#!/usr/bin/env bash
# =============================================================================
# Round-7 MC trigger-efficiency regeneration driver
# (docs/tracking/mc_trigger_efficiency.md, round-7 Autonomy Contract item 1).
#
# Regenerates the COMPLETE MC trig-eff product set on top of the two round-7 code changes:
#   (a) the TRUTH FIDUCIAL (truth pT > 4.5, |truth eta| < 2.4) is now part of the sample
#       selection for every sample and every step (FillMCTrigEffHists.cxx);
#   (b) the Step-3/Step-4 ratio ERROR BARS use the conditional (binomial-correct) form
#       instead of TH1::Divide's independent propagation (errA/errB/covP/covQ histograms
#       + SetConditionalRatioErrors in plot_mc_trig_eff.cxx). Central values unchanged.
#
# Matrix: {pp_full, overlay, noovl} x {Tight, Medium} x {step1/2, fit, step3, step4} -> plots.
#
# CONCURRENCY: the three ACLiC macros are pre-compiled ONCE (a race on the shared _cxx.so is
# the one real hazard), after which the per-sample fill/fit chains are independent -- they
# read different input files and write different output files -- so they run in parallel.
# The PLOTTING is serialized at the end because "noovl" makes a three-way comparison and
# needs the pp_full and overlay fits to exist.
#
# Usage:  ./run_mc_trigeff_round7.sh            (everything)
#         SAMPLES="overlay noovl" ./run_mc_trigeff_round7.sh
#         WPS="tight" ./run_mc_trigeff_round7.sh
# =============================================================================
set -Euo pipefail

SAMPLES="${SAMPLES:-pp_full overlay noovl}"
WPS="${WPS:-tight medium}"

ANALYSIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
PLOT_DIR="${ANALYSIS_DIR}/plotting_codes/trig_effcy/mc_based"
LOG_DIR="${ANALYSIS_DIR}/pipelines/logs_round7"
mkdir -p "$LOG_DIR"

ts(){ date "+%Y-%m-%d %H:%M:%S"; }
log(){ echo "[$(ts)] $*"; }
fail(){ echo "[$(ts)] FATAL: $*" >&2; exit 1; }

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
set +eu
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh" --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"
set -eu

# sample -> output directory + hist-file label (mirrors FillMCTrigEffHists::GetSampleConfig)
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

# A ROOT macro that throws still exits 0 -- never trust the exit code, check the artefact.
val_file(){ [[ -s "$1" ]] || fail "missing/empty artefact: $1"; }

# ---- 0. pre-compile the three macros ONCE (serial: shared _cxx.so) ----------------------
log "════ pre-compiling ACLiC macros ════"
# NB: `root -q 'X.cxx+'` compiles AND RUNS X() with default arguments -- use `.L X.cxx+`
# so this stage only builds the shared library.
( cd "$RDF_DIR"  && root -l -b -q -e '.L FillMCTrigEffHists.cxx+' ) >"$LOG_DIR/compile_fill.log" 2>&1 \
    || fail "FillMCTrigEffHists compile failed (see $LOG_DIR/compile_fill.log)"
( cd "$RDF_DIR"  && root -l -b -q -e '.L FitMCSinglesEffcy.cxx+'  ) >"$LOG_DIR/compile_fit.log"  2>&1 \
    || fail "FitMCSinglesEffcy compile failed"
( cd "$PLOT_DIR" && root -l -b -q -e '.L plot_mc_trig_eff.cxx+'   ) >"$LOG_DIR/compile_plot.log" 2>&1 \
    || fail "plot_mc_trig_eff compile failed"
( cd "$PLOT_DIR" && root -l -b -q -e '.L plot_mc_singles_2d_effcy.cxx+' ) \
    >"$LOG_DIR/compile_singles2d.log" 2>&1 || fail "plot_mc_singles_2d_effcy compile failed"
( cd "$PLOT_DIR" && root -l -b -q -e '.L write_mc_pair_statistics_tables.cxx+' ) \
    >"$LOG_DIR/compile_stats.log" 2>&1 || fail "write_mc_pair_statistics_tables compile failed"
log "  compiled OK"

# ---- 1. fill / fit chain, one background job per (sample, WP) ---------------------------
chain(){   # $1=sample $2=wp
    local s="$1" wp="$2" cpp suf out lbl lg
    cpp=$([[ $wp == tight ]] && echo true || echo false)
    suf=$([[ $wp == tight ]] && echo ""   || echo "_medium_wp")
    out="$(outdir_of "$s")"; lbl="$(label_of "$s")"
    lg="$LOG_DIR/${s}_${wp}.log"
    {
        echo "=== $s / $wp : step1+2 ==="
        ( cd "$RDF_DIR" && root -l -b -q "FillMCTrigEffHists.cxx+(\"$s\", false, $cpp)" )
        val_file "$out/mc_trig_eff_hists_${lbl}${suf}.root"
        echo "=== $s / $wp : fit ==="
        ( cd "$RDF_DIR" && root -l -b -q "FitMCSinglesEffcy.cxx+(\"$s\", $cpp)" )
        val_file "$out/single_mu_effcy_pT_fit_mc${suf}.root"
        echo "=== $s / $wp : step3 ==="
        ( cd "$RDF_DIR" && root -l -b -q "FillMCTrigEffHists.cxx+(\"$s\", true, $cpp)" )
        val_file "$out/mc_trig_eff_hists_${lbl}${suf}_step3.root"
        echo "=== $s / $wp : step4 ==="
        ( cd "$RDF_DIR" && root -l -b -q "FillMCTrigEffHists.cxx+(\"$s\", false, $cpp, true)" )
        val_file "$out/mc_trig_eff_hists_${lbl}${suf}_step4.root"
        # Step-1 SANITY CHECK (Physics Procedure 3.5). Without this stage a clean regeneration
        # would silently produce no sanity plots at all: the plot macro skips the block with a
        # note when the _sanity.root is absent, so its absence is NOT an error anywhere else.
        echo "=== $s / $wp : sanity ==="
        ( cd "$RDF_DIR" && root -l -b -q "FillMCTrigEffHists.cxx+(\"$s\", false, $cpp, false, true)" )
        val_file "$out/mc_trig_eff_hists_${lbl}${suf}_sanity.root"
        echo "=== $s / $wp : CHAIN OK ==="
    } >"$lg" 2>&1
}

pids=(); names=()
for s in $SAMPLES; do for wp in $WPS; do
    log "launching chain: $s / $wp  -> $LOG_DIR/${s}_${wp}.log"
    chain "$s" "$wp" & pids+=($!); names+=("$s/$wp")
done; done

rc=0
for i in "${!pids[@]}"; do
    if wait "${pids[$i]}"; then log "  chain OK : ${names[$i]}"
    else log "  chain FAILED : ${names[$i]}"; rc=1; fi
done
[[ $rc -eq 0 ]] || fail "at least one fill/fit chain failed -- see $LOG_DIR/*.log"

# ---- 2. plots (serial; noovl is a three-way comparison and needs the other two fits) ----
for wp in $WPS; do
    cpp=$([[ $wp == tight ]] && echo true || echo false)
    for s in $SAMPLES; do
        log "plotting $s / $wp"
        ( cd "$PLOT_DIR" && root -l -b -q "plot_mc_trig_eff.cxx+(\"$s\", $cpp)" ) \
            >"$LOG_DIR/plot_${s}_${wp}.log" 2>&1 || fail "plot $s/$wp failed"
        grep -q "done\." "$LOG_DIR/plot_${s}_${wp}.log" || fail "plot $s/$wp did not reach the end"
        grep -q "SANITY CHECK" "$LOG_DIR/plot_${s}_${wp}.log" \
            || fail "plot $s/$wp produced no Step-1 sanity block (missing _sanity.root?)"

        # Step-1 2D efficiency maps (round 9). r17663 is a one-time Step-1 cross-check and
        # deliberately gets no analysis-level plot set, so it is skipped here.
        #
        # The exit code proves nothing (see the val_file comment above): a throwing macro still
        # exits 0. The ARTEFACT is the pair of PNGs, and their directory is built in C++ from
        # DrCorrSampleCfg::out_base, which this shell must NOT re-derive (the dr-correction driver
        # records what happens when the two constructions drift). So take the paths from the
        # macro's own two "wrote" lines and check the files they name.
        if [[ $s != noovl ]]; then
            s2d_log="$LOG_DIR/singles2d_${s}_${wp}.log"
            ( cd "$PLOT_DIR" && root -l -b -q "plot_mc_singles_2d_effcy.cxx+(\"$s\", $cpp)" ) \
                >"$s2d_log" 2>&1 || fail "2D singles plot $s/$wp failed"
            mapfile -t s2d_png < <(sed -n 's/^ *wrote //p' "$s2d_log")
            [[ ${#s2d_png[@]} -eq 2 ]] \
                || fail "2D singles plot $s/$wp wrote ${#s2d_png[@]} PNG(s), expected 2 -- see $s2d_log"
            for f in "${s2d_png[@]}"; do val_file "$f"; done
            for want in step1_eff_2d_pt_vs_q_eta_charge_sepr.png \
                        step1_eff_2d_pt_vs_q_eta_charge_comb.png; do
                printf '%s\n' "${s2d_png[@]}" | grep -q "/step1_singles_data_mc/${want}\$" \
                    || fail "2D singles plot $s/$wp did not write $want -- see $s2d_log"
            done
        fi

        # Muon-pair count / cross-section tables (round 9). ONE failure mode is benign and
        # expected during the transition: a sample filled BEFORE the round-9 STATISTICS
        # BOOKKEEPING block simply has no statistics histograms, and the macro says so with
        # "histogram '...' is MISSING from <file>". Every OTHER throw is a hard guard --
        # pair-pT bin count vs MCTrigEffPairPt::NBins(), axis-edge mismatch among the six TH2Ds,
        # an all-empty binned range -- i.e. a binning or provenance error, and downgrading those
        # to a log line is exactly how a wrong table gets published. So: note only the MISSING
        # case, fail on anything else.
        if [[ $s != noovl ]]; then
            st_log="$LOG_DIR/stats_tables_${s}_${wp}.log"
            ( cd "$PLOT_DIR" && root -l -b -q "write_mc_pair_statistics_tables.cxx+(\"$s\", $cpp)" ) \
                >"$st_log" 2>&1 || true
            if grep -q "CSVs written to" "$st_log"; then
                :
            elif grep -q "is MISSING from" "$st_log"; then
                log "  note: no statistics tables for $s/$wp (needs a post-round-9 Step-3 fill)"
            else
                fail "statistics tables $s/$wp failed for a reason other than a pre-round-9 fill -- see $st_log"
            fi
        fi
    done
done

log "════ round-7 regeneration DONE ════"
