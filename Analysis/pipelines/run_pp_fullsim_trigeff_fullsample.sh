#!/usr/bin/env bash
set -Euo pipefail
# Task-2 re-run driver: pp24 FULL-sample MC trigger efficiency ONLY.
# The full production run already did the nominal + single-muon NTP and all non-trig-eff plots;
# its Stage 3 (mc_trig NTP) hit the store_mc_trigger multi-file bug (fixed in 31c2724). This
# redoes just the two mc_trig NTP passes + the Stage-10 chain, so nominal/single NTP is NOT redone.
#
# Hardening lesson from that bug: a ROOT stage that throws still exits 0, so `|| fail` misses it.
# Here every NTP/step output is VALIDATED to exist and be non-empty before proceeding.
#
# Usage: USE_TIGHT_WP=1 ./run_pp_fullsim_trigeff_fullsample.sh

USE_TIGHT_WP="${USE_TIGHT_WP:-1}"
TIGHT_CPP=$([[ $USE_TIGHT_WP == 1 ]] && echo true || echo false)

ANALYSIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NTP_DIR="${ANALYSIS_DIR}/NTupleProcessingCode"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
PLOT_DIR="${ANALYSIS_DIR}/plotting_codes/trig_effcy/mc_based"
FULL=/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample
PP_TRIG=/usatlas/u/yuhanguo/usatlasdata/dimuon_data/plots/pp_trigger_efficiency/mc_based

MCTRIG_PAIR="${FULL}/muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_full.root"
MCTRIG_SINGLE="${FULL}/muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_mc_trig_single_muon_full.root"

ts(){ date "+%Y-%m-%d %H:%M:%S"; }
log(){ echo "[$(ts)] $*"; }
fail(){ echo "[$(ts)] FATAL: $*" >&2; exit 1; }

# validate a tree is present and non-empty (ROOT exits 0 on a swallowed throw, so a script's
# exit code is NOT enough -- check the artefact).
val_tree(){  # $1=file $2=tree
    [[ -f "$1" ]] || fail "missing $1"
    local n
    n=$(root -l -b -q -e "TFile*f=TFile::Open(\"$1\"); TTree*t=f&&!f->IsZombie()?(TTree*)f->Get(\"$2\"):nullptr; printf(\"N %lld\n\", t?t->GetEntries():0LL);" 2>/dev/null | grep -oP 'N \K[0-9]+' | tail -1 || true)
    [[ "${n:-0}" -gt 0 ]] || fail "$1 : tree $2 is empty/missing (store_mc_trigger silent throw?)"
    log "  validated $(basename "$1") : $2 = ${n} entries"
}

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
set +eu
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh" --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"
set -eu

log "════ pp24 FULL-sample MC trig-eff re-run (WP=$([[ $USE_TIGHT_WP == 1 ]] && echo tight || echo medium)) ════"

# --- mc_trig NTP (two full passes over the farm) ---
log "[1] mc_trig pair NTP"
( cd "$NTP_DIR" && bash run_pythia_fullsim_mc_trig_full_sample.sh ) || fail "mc_trig pair NTP script failed"
val_tree "$MCTRIG_PAIR" "muon_pair_tree_kin0_sign2"

log "[2] mc_trig single-muon NTP"
( cd "$NTP_DIR" && bash run_pythia_fullsim_single_muon_mc_trig_full_sample.sh ) || fail "mc_trig single-muon NTP script failed"
val_tree "$MCTRIG_SINGLE" "muon_tree"

# --- Stage 10 chain (pp_full) ---
if [[ -d "$PP_TRIG" ]]; then
    bak="${PP_TRIG}.bak_$(date +%Y%m%d_%H%M%S)"; cp -a "$PP_TRIG" "$bak"; log "  backed up canonical pp trig-eff -> $(basename "$bak")"
fi

log "[3] FillMCTrigEffHists step1 (pp_full)"
( cd "$RDF_DIR" && root -l -b -q 'FillMCTrigEffHists.cxx+("pp_full", false, '"$TIGHT_CPP"')' ) || fail "FillMCTrigEffHists step1 failed"
HISTS="${FULL}/mc_trig_eff_hists_pp24_full$([[ $USE_TIGHT_WP == 1 ]] || echo _medium_wp).root"
[[ -f "$HISTS" ]] || fail "step1 produced no hist file $HISTS"
log "  step1 hists: $(basename "$HISTS")"

log "[4] FitMCSinglesEffcy (pp_full)"
( cd "$RDF_DIR" && root -l -b -q 'FitMCSinglesEffcy.cxx+("pp_full", '"$TIGHT_CPP"')' ) || fail "FitMCSinglesEffcy failed"

log "[5] FillMCTrigEffHists step3 (pp_full)"
( cd "$RDF_DIR" && root -l -b -q 'FillMCTrigEffHists.cxx+("pp_full", true, '"$TIGHT_CPP"')' ) || fail "FillMCTrigEffHists step3 failed"

log "[6] plot_mc_trig_eff (pp_full)"
( cd "$PLOT_DIR" && root -l -b -q 'plot_mc_trig_eff.cxx+("pp_full", '"$TIGHT_CPP"')' ) || fail "plot_mc_trig_eff failed"

log "════ DONE -> plots in ${PP_TRIG} ════"
