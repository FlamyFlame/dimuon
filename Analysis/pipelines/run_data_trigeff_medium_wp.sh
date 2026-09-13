#!/usr/bin/env bash
set -Euo pipefail
# =============================================================================================
# MEDIUM-WP pass of the DATA single-muon mu4 trigger efficiency (tag-and-probe), pp24 + PbPb
# 23/24/25/26.
#
# WHY THIS SCRIPT EXISTS. `pipeline_pp_trig_eff.sh` and `pipeline_pbpb_trig_eff.sh` run the
# NOMINAL working point only (`isTight` defaults to true), so every rerun of those pipelines
# refreshes the Tight tag-and-probe file and silently leaves the Medium one behind. That matters
# because the Medium file is not a spare: `plot_mc_trig_eff.cxx` with `use_tight_wp = false`
# reads the `_medium_wp` DATA file as the reference the Medium MC is compared against (it never
# compares across working points -- see its header). A stale Medium file therefore produces a
# comparison in which the data and the MC had DIFFERENT selections applied, with no error and no
# warning. The repo requires every plot set and its code to expose a Tight/Medium switch
# (`Analysis/docs/muon_wp_registry.md`), and this is the data trig-eff half of it.
#
# WHAT IT RUNS, per sample: exactly the two WP-dependent stages of Pipeline 2 --
#   (1) RDF tag-and-probe hist filling with isTight = false  -> ..._qeta_fid_medium_wp.root
#   (2) the pT turn-on fit with wp_suffix = "_medium_wp"     -> single_mu_effcy_pT_fit_medium_wp.root
# and nothing else. There is no Medium DATA plot set (the data trig-eff plots are Tight-only),
# so no plotting stage: the Medium products exist to be consumed by the Medium MC comparison.
#
# The NTuple processing is NOT rerun and does not need to be: the working point and the q*eta
# fiducial gap cut are both applied at the RDF stage, on branches the NTuple processing already
# stores (`pair_pass_tight` / `pair_pass_medium`, `q_eta2nd`).
#
# Usage:   ./run_data_trigeff_medium_wp.sh
#          SAMPLES="pp 23 24 25 26" ./run_data_trigeff_medium_wp.sh
#          RDF_NTHREADS=4 ./run_data_trigeff_medium_wp.sh
# =============================================================================================

# Pb+Pb years here are FILTERED on data presence, so the default can safely name every
# period: a year whose raw skim NTUPs are not on disk yet is announced and skipped rather
# than failing the run after all the real work is done.  "pp" is never filtered.
SAMPLES="${SAMPLES:-pp 23 24 25 26}"
_data_base="/usatlas/u/yuhanguo/usatlasdata/dimuon_data"
_kept=()
for _s in ${SAMPLES}; do
  # Strict match -- see the note in pipeline_pbpb_crossx.sh (.bak re-merge files).
  if [[ "${_s}" == "pp" ]] || ls "${_data_base}/pbpb_20${_s}" 2>/dev/null \
       | grep -qE "^data_pbpb${_s}_part[0-9]+\.root$"; then
    _kept+=("${_s}")
  else
    echo "[SKIP] Pb+Pb 20${_s}: no raw skim NTUPs on disk -- sample skipped." >&2
  fi
done
if [[ ${#_kept[@]} -eq 0 ]]; then
  echo "[FATAL] no requested sample has input on disk." >&2
  exit 1
fi
SAMPLES="${_kept[*]}"
echo "[INFO] samples to process: ${SAMPLES}"
RDF_NTHREADS="${RDF_NTHREADS:-2}"

ANALYSIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RDF_DIR="${ANALYSIS_DIR}/RDFBasedHistFilling"
DATA_BASE=/usatlas/u/yuhanguo/usatlasdata/dimuon_data

ts(){ date "+%Y-%m-%d %H:%M:%S"; }
log(){ echo "[$(ts)] $*"; }
fail(){ echo "[$(ts)] FATAL: $*" >&2; exit 1; }

# A ROOT macro that throws still exits 0 (the exception aborts the interpreter after ROOT has
# decided the job "ran"), so every stage below is judged by the artefact it must have written,
# never by the exit code -- the same rule the two pipelines this complements use.
val_new(){   # $1 = file that must exist, be non-empty, and be NEWER than the marker file $2
  [[ -s "$1" ]] || fail "missing/empty artefact: $1"
  [[ "$1" -nt "$2" ]] || fail "artefact not rewritten by this run (older than the run marker): $1"
  log "  OK $(basename "$1")"
}

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
set +eu
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh" --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"
set -eu

MARKER="$(mktemp)"; trap 'rm -f "${MARKER}"' EXIT

log "════ DATA single-muon mu4 trigger efficiency, MEDIUM WP: ${SAMPLES} ════"

for s in ${SAMPLES}; do
  if [[ "$s" == "pp" ]]; then
    log "[pp24] RDF tag-and-probe fill (Medium WP)"
    ( cd "$RDF_DIR" && root -l -b <<EOF
ROOT::EnableImplicitMT(${RDF_NTHREADS});
.L RDFBasedHistFillingPP.cxx+
RDFBasedHistFillingPP pp(24);
pp.trigger_mode = 1;
pp.isTight = false;          // the ONLY difference from the nominal pipeline stage
pp.Run();
gSystem->Exit(0);
EOF
    ) || fail "pp24 Medium RDF fill returned non-zero"
    val_new "${DATA_BASE}/pp_2024/histograms_real_pairs_pp_2024_single_mu4_coarse_q_eta_bin_qeta_fid_medium_wp.root" "${MARKER}"

    log "[pp24] pT turn-on fit (erf+log, Medium WP)"
    ( cd "$ANALYSIS_DIR" && root -l -b <<EOF
.L SingleMuEffcyPtTurnOnFitter.cxx
single_muon_trig_effcy_pT_fitting("_medium_wp");
gSystem->Exit(0);
EOF
    ) || fail "pp24 Medium turn-on fit returned non-zero"
    val_new "${DATA_BASE}/pp_2024/trg_effcy_pT_fitting_to_erf_plus_log/single_mu_effcy_pT_fit_medium_wp.root" "${MARKER}"
  else
    yr="$s"
    log "[PbPb 20${yr}] RDF tag-and-probe fill (Medium WP)"
    ( cd "$RDF_DIR" && root -l -b <<EOF
ROOT::EnableImplicitMT(${RDF_NTHREADS});
.L RDFBasedHistFillingPbPb.cxx+
RDFBasedHistFillingPbPb pbpb(${yr});
pbpb.trigger_mode = 1;
pbpb.mindR_trig = 0.02;
pbpb.isTight = false;        // the ONLY difference from the nominal pipeline stage
pbpb.Run();
gSystem->Exit(0);
EOF
    ) || fail "PbPb 20${yr} Medium RDF fill returned non-zero"
    val_new "${DATA_BASE}/pbpb_20${yr}/histograms_real_pairs_pbpb_20${yr}_single_mu4_coarse_q_eta_bin_qeta_fid_medium_wp.root" "${MARKER}"

    log "[PbPb 20${yr}] pT turn-on fit (Fermi+log, Medium WP)"
    ( cd "$ANALYSIS_DIR" && root -l -b <<EOF
.L SingleMuEffcyPtTurnOnFitter.cxx
single_muon_trig_effcy_pT_fitting_PbPb(${yr}, "_medium_wp");
gSystem->Exit(0);
EOF
    ) || fail "PbPb 20${yr} Medium turn-on fit returned non-zero"
    val_new "${DATA_BASE}/pbpb_20${yr}/trg_effcy_pT_fitting_to_fermi_plus_log/single_mu_effcy_pT_fit_medium_wp.root" "${MARKER}"
  fi
done

log "════ Medium-WP data trigger efficiency DONE ════"
