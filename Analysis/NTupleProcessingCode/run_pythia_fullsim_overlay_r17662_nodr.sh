#!/bin/bash
# r17662 (signal-only-truth) overlay NTP, pTH8_14, FULL 10k events.
# Runs the standard procedure: pure prob>0.5 barcode truth->reco matching (the ad-hoc
# dR<0.05 fallback was deleted 2026-07-10 — it over-estimated the reco efficiency for
# ΔR<0.05 pairs; see hijing_overlay_det_response_band.md).
#
# Input: r17662_run/ contains a symlink of the r17662 NTUP under the expected
#   FullSimHIJINGOverlayPP24 name; the other 5 pT slices are absent (skipped).
#
# Usage: ./run_pythia_fullsim_overlay_r17662_nodr.sh [single_muon|pair]
MODE="${1:-single_muon}"
cd "$(dirname "$0")"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

INDIR=/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/r17662_run/

if [[ "$MODE" == "single_muon" ]]; then
  SINGLE="true"; SUFFIX="_r17662_nodr_single_muon"
else
  SINGLE="false"; SUFFIX="_r17662_nodr"
fi

root -b -l << EOF
	.L PythiaAnalysisClasses.h+

	PythiaFullSimOverlayAnalysis py(FullSimSampleType::hijing);
	py.fullsim_input_dir_override = "${INDIR}";
	py.fill_kn_trees_fullsim = true;
	py.output_single_muon_tree = ${SINGLE};
	py.extra_output_suffix = "${SUFFIX}";
	py.Run();

	.q;

EOF
