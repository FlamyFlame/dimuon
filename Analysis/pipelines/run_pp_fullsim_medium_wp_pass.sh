#!/usr/bin/env bash
set -Euo pipefail
# MEDIUM-WP systematic pass for the pp24 FULL sample. NTP + RDF hists are WP-agnostic (they store
# both pass_medium and pass_tight), so this only re-runs the PLOTTERS/FITS at Medium WP.
#
# Runs ONLY the plots whose output is WP-separated (so Tight nominal is never clobbered):
#   - reco efficiency + detector response  -> pp24_reco_effcy_plots/medium/ , pp24_det_resp_plots/medium/
#   - single-muon reco efficiency          -> pp24_single_muon_reco_effcy/medium_wp/
#   - MC trigger efficiency                -> pp_trigger_efficiency/mc_based/step*/medium/
# DELIBERATELY EXCLUDED: crossx / kn table / reco-distr -- those are nominal cross-section /
# statistics deliverables with a single (Tight) WP and NO tight/medium output separation; running
# them at Medium would overwrite the validated Tight nominal.
#
# Prereq: the FULL-sample NTP (nominal + single + mc_trig) and RDF hists already exist (task 1/2
# Tight run). Usage: ./run_pp_fullsim_medium_wp_pass.sh

ANALYSIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RECO_DIR="${ANALYSIS_DIR}/plotting_codes/reco_effcy"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
TRIG_PLOT_DIR="${ANALYSIS_DIR}/plotting_codes/trig_effcy/mc_based"
PP_TRIG=/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/mc_based

ts(){ date "+%Y-%m-%d %H:%M:%S"; }
log(){ echo "[$(ts)] $*"; }
fail(){ echo "[$(ts)] FATAL: $*" >&2; exit 1; }

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
set +eu
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh" --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"
set -eu

log "════ pp24 FULL-sample MEDIUM-WP pass ════"

log "[1] reco efficiency + detector response (medium)"
( cd "$RECO_DIR" && root -l -b <<ROOTEOF
.L PythiaFullsimRecoEffPlotter.cxx+
{ PythiaFullsimRecoEffPlotter pl(false); pl.is_test_sample = false; pl.Run(); }
.q
ROOTEOF
) || fail "reco-eff/det-resp medium failed"

log "[2] single-muon reco efficiency (medium)"
( cd "$RECO_DIR" && root -l -b <<ROOTEOF
.L plot_single_muon_reco_effcy.cxx+
plot_single_muon_reco_effcy("pp", "default", false, false);
.q
ROOTEOF
) || fail "single-muon reco-eff medium failed"

# --- MC trigger efficiency, medium ---
if [[ -d "$PP_TRIG" ]]; then
    bak="${PP_TRIG}.bak_$(date +%Y%m%d_%H%M%S)"; cp -a "$PP_TRIG" "$bak"; log "  backed up canonical pp trig-eff -> $(basename "$bak")"
fi
log "[3] FillMCTrigEffHists step1 (pp_full, medium)"
( cd "$RDF_DIR" && root -l -b -q 'FillMCTrigEffHists.cxx+("pp_full", false, false)' ) || fail "Fill step1 medium failed"
log "[4] FitMCSinglesEffcy (pp_full, medium)"
( cd "$RDF_DIR" && root -l -b -q 'FitMCSinglesEffcy.cxx+("pp_full", false)' ) || fail "Fit medium failed"
log "[5] FillMCTrigEffHists step3 (pp_full, medium)"
( cd "$RDF_DIR" && root -l -b -q 'FillMCTrigEffHists.cxx+("pp_full", true, false)' ) || fail "Fill step3 medium failed"
log "[6] plot_mc_trig_eff (pp_full, medium)"
( cd "$TRIG_PLOT_DIR" && root -l -b -q 'plot_mc_trig_eff.cxx+("pp_full", false)' ) || fail "plot_mc_trig_eff medium failed"

log "════ MEDIUM-WP pass DONE ════"
