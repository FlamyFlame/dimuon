#!/bin/bash
# pbpb23 (r17618) overlay NTP restricted to the pT-hat 125-300 GeV slice -- the like-for-like
# reference for the single-slice Pb+Pb-2024-conditions test sample (r17864, pTH125_300 only;
# docs/tracking/hijing_overlay_pbpb24_test_sample_skim.md, verification round). Same NTP code,
# same date, same slice, so only the r-tag differs.
#
# Input: r17618_pTH125_300_run/ holds a symlink to the pTH125_300 NTUP of the _pbpb23 sample
#   under its expected name; the other 5 slices are absent (skipped). Outputs land in that
#   directory (fullsim_input_dir_override => output_dir), never on the nominal 6-slice products.
#
# Usage: ./run_pythia_fullsim_overlay_pbpb23_pTH125_300_slice.sh [pair|single_muon|mc_trig|mc_trig_single_muon]
MODE="${1:-pair}"
cd "$(dirname "$0")"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"
# pinned to the _pbpb23 directory explicitly (FullSimSampleType.h: the unsuffixed dir is pbpb24)
INDIR=/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample_pbpb23/r17618_pTH125_300_run/
case "$MODE" in
  pair)                SINGLE=false; TRIG=false; SUFFIX="_r17618pTH125_300" ;;
  single_muon)         SINGLE=true;  TRIG=false; SUFFIX="_r17618pTH125_300_single_muon" ;;
  mc_trig)             SINGLE=false; TRIG=true;  SUFFIX="_r17618pTH125_300" ;;              # store_mc_trigger adds _mc_trig itself
  mc_trig_single_muon) SINGLE=true;  TRIG=true;  SUFFIX="_r17618pTH125_300_single_muon" ;;  # -> ..._mc_trig_r17618pTH125_300_single_muon
  *) echo "unknown mode $MODE"; exit 1 ;;
esac
root -b -l << EOF
	.L PythiaAnalysisClasses.h+
	PythiaFullSimOverlayAnalysis py(FullSimSampleType::hijing);
	// TEST sample: reads the overlay test sample AND uses the pp beam only (only beam produced).
	// ONE switch -- it also selects the input dir. See FullSimSampleType.h.
	py.isTestSample = true;
	py.overlay_pbpb_year = 23;   // label hijing_overlay_pbpb23, matching INDIR above
	// AMI PROVENANCE: the overlay is built on the pp-beam evgen DSIDs. AMI files are keyed by
	// beam+slice only, so declare them; InitInputFullsim throws if another production is read.
	// (Overlay -> nPDF evgen -> the e8599 files in the evgen's ami_info_nPDF/, FullSimSampleType.h.)
	py.expected_ami_dsids = {802776, 802777, 802778, 802779, 802780, 802781};
	// Single-slice DIAGNOSTIC run: INDIR holds pTH125_300 ONLY. Every output is a ratio or a
	// single-slice shape, so a missing slice is intentional and harmless here.
	py.allow_missing_slices = true;
	py.fullsim_input_dir_override = "${INDIR}";
	py.fill_kn_trees_fullsim = true;
	py.store_mc_trigger = ${TRIG};
	py.output_single_muon_tree = ${SINGLE};
	py.extra_output_suffix = "${SUFFIX}";
	py.Run();
	.q;
EOF
