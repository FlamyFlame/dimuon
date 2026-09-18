#!/bin/bash
# fullsim_sample_layout.sh -- SHELL TWIN of MuonObjectsParamsAndHelpers/FullSimSampleType.h
# ("PER-SAMPLE DIRECTORY LAYOUT") and of dr_correction_sample_cfg.h (sample keys -> dir / label).
#
# The C++ header is the source of truth; this file exists only because bash pipelines have to
# LOCATE the products the macros write (freshness checks, validation, backups) and cannot include
# the header. Every subdirectory name and every basename prefix below is a verbatim copy of the
# header's -- when one changes, change both in the same commit. `source` this file; never copy
# its strings into another script.
#
#   <sample>/                      raw NTUP, merging-record.txt, grid-monitor state (NO ami_info/:
#                                  AMI weights are keyed by the Pythia EVGEN, FullSimSampleType.h)
#     muon_pairs_*.root, hists_pythia_ntuple_processing_*.root, histograms_pythia_fullsim_*.root
#                                  (FLAT, like the data directories)
#     mc_trig_eff/hists/           mc_trig_eff_hists_<label>*.root
#     mc_trig_eff/singles_fits/    single_mu_effcy_pT_fit_mc_<label>*.root
#     mc_trig_eff/dr_correction/   dr_correction_plateaus_* , dr_correction_fits_*
#     mc_trig_eff/pair_eff/        pair_trig_eff_<label>*.root
#     mc_trig_eff/closure/         mc_trig_eff_closure_<label>*.root
#     reco_eff/                    pair_reco_eff_<label>.root
#     plots/  backup/  logs/
#
# Sample keys (dr_correction_sample_cfg.h::GetDrCorrSample): pp | pp_full | overlay | noovl.
# OVERLAY_YEAR (env, default 24) = the HIJING-overlay conditions year the "overlay" key means:
#   24 -> pythia_fullsim_hijing_overlay_test_sample/         label hijing_overlay_pbpb24
#   23 -> pythia_fullsim_hijing_overlay_test_sample_pbpb23/  label hijing_overlay_pbpb23

FULLSIM_DATA_ROOT="/usatlas/u/yuhanguo/usatlasdata"
: "${OVERLAY_YEAR:=24}"

fullsim_check_overlay_year() {
  case "${1:-$OVERLAY_YEAR}" in
    23|24) ;;
    *) echo "fullsim_sample_layout.sh: OVERLAY_YEAR must be 23 or 24, got '${1:-$OVERLAY_YEAR}'" >&2; return 1 ;;
  esac
}

# sample key -> sample ROOT directory (trailing slash), mirrors FullSimSampleInputDir
fullsim_sample_dir() {
  local key="$1" year="${2:-$OVERLAY_YEAR}"
  case "$key" in
    pp)      echo "${FULLSIM_DATA_ROOT}/pythia_fullsim_test_sample/" ;;
    pp_full) echo "${FULLSIM_DATA_ROOT}/pythia_fullsim_full_sample/" ;;
    overlay)
      fullsim_check_overlay_year "$year" || return 1
      if [[ "$year" == "24" ]]; then echo "${FULLSIM_DATA_ROOT}/pythia_fullsim_hijing_overlay_test_sample/"
      else echo "${FULLSIM_DATA_ROOT}/pythia_fullsim_hijing_overlay_test_sample_pbpb23/"; fi ;;
    noovl)   echo "${FULLSIM_DATA_ROOT}/pythia_fullsim_no_overlay_test_sample/" ;;
    *) echo "fullsim_sample_layout.sh: unknown sample key '$key'" >&2; return 1 ;;
  esac
}

# sample key -> file label, mirrors GetDrCorrSample / FullSimSampleLabel
fullsim_sample_label() {
  local key="$1" year="${2:-$OVERLAY_YEAR}"
  case "$key" in
    pp)      echo "pp24" ;;
    pp_full) echo "pp24_full" ;;
    overlay) fullsim_check_overlay_year "$year" || return 1; echo "hijing_overlay_pbpb${year}" ;;
    noovl)   echo "r17663_no_overlay" ;;
    *) echo "fullsim_sample_layout.sh: unknown sample key '$key'" >&2; return 1 ;;
  esac
}

# subtrees of a sample directory (argument = sample dir WITH trailing slash)
fullsim_mc_trig_eff_hists_dir()   { echo "${1}mc_trig_eff/hists/"; }
fullsim_mc_trig_eff_fits_dir()    { echo "${1}mc_trig_eff/singles_fits/"; }
fullsim_mc_trig_eff_drcorr_dir()  { echo "${1}mc_trig_eff/dr_correction/"; }
fullsim_mc_trig_eff_paireff_dir() { echo "${1}mc_trig_eff/pair_eff/"; }
fullsim_mc_trig_eff_closure_dir() { echo "${1}mc_trig_eff/closure/"; }
fullsim_reco_eff_dir()            { echo "${1}reco_eff/"; }
fullsim_plots_dir()               { echo "${1}plots/"; }
fullsim_backup_dir()              { echo "${1}backup/"; }
fullsim_logs_dir()                { echo "${1}logs/"; }

# single_mu_effcy_pT_fit_mc_<label><corr_suf><wp_suf>.root  (FullSimMCSinglesFitFile)
#   $1 sample dir, $2 label, $3 corr suffix ("" | _corrected | _corrected_sfclosure), $4 wp suffix ("" | _medium_wp)
fullsim_singles_fit_file() {
  echo "$(fullsim_mc_trig_eff_fits_dir "$1")single_mu_effcy_pT_fit_mc_${2}${3:-}${4:-}.root"
}
# mc_trig_eff_hists_<label><wp_suf><variant>.root  (DrCorrHistFile)
fullsim_hists_file() {
  echo "$(fullsim_mc_trig_eff_hists_dir "$1")mc_trig_eff_hists_${2}${3:-}${4:-}.root"
}
# mc_trig_eff_closure_<label><wp_suf><variant>.root  (DrCorrClosureFile)
fullsim_closure_file() {
  echo "$(fullsim_mc_trig_eff_closure_dir "$1")mc_trig_eff_closure_${2}${3:-}${4:-}.root"
}
# pair_trig_eff_<label><wp_suf>.root  (PairTrigEff::FileName)
fullsim_pair_trig_eff_file() {
  echo "$(fullsim_mc_trig_eff_paireff_dir "$1")pair_trig_eff_${2}${3:-}.root"
}
# dr_correction_plateaus_<label><wp_suf><ptbin_suf>.root  (DrCorrPlateauFile)
fullsim_plateau_file() {
  echo "$(fullsim_mc_trig_eff_drcorr_dir "$1")dr_correction_plateaus_${2}${3:-}${4:-}.root"
}
# dr_correction_fits_<label><wp_suf><ptbin_suf>_step<N>_<method><sign_suf><mode_tag>.root  (DrCorrFitFile)
#   $1 dir $2 label $3 wp_suf $4 ptbin_suf $5 step $6 method $7 sign_suf $8 mode_tag
fullsim_fit_file() {
  echo "$(fullsim_mc_trig_eff_drcorr_dir "$1")dr_correction_fits_${2}${3:-}${4:-}_step${5}_${6}${7:-}${8:-}.root"
}

# Back up a product before it is overwritten: <sample>/backup/<basename>.bak_<stamp>.root
#   $1 file, $2 sample dir. Prints the backup path.
fullsim_backup_file() {
  local f="$1" sdir="$2" bdir stamp
  bdir="$(fullsim_backup_dir "$sdir")"; mkdir -p "$bdir"
  stamp="$(date +%Y%m%d_%H%M%S)"
  local base; base="$(basename "$f")"
  local bak="${bdir}${base%.root}.bak_${stamp}.root"
  cp -a "$f" "$bak" && echo "$bak"
}
