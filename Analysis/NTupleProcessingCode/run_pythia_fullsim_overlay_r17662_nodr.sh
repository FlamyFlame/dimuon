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

# r17662 is a Pb+Pb-2023-conditions reprocessing of the pbpb23 test sample -> pinned to the
# _pbpb23 directory explicitly (FullSimSampleType.h: the unsuffixed dir is the pbpb24 sample).
INDIR=/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_hijing_overlay_test_sample_pbpb23/r17662_run/

if [[ "$MODE" == "single_muon" ]]; then
  SINGLE="true"; SUFFIX="_r17662_nodr_single_muon"
else
  SINGLE="false"; SUFFIX="_r17662_nodr"
fi

root -b -l << EOF
	.L PythiaAnalysisClasses.h+

	PythiaFullSimOverlayAnalysis py(FullSimSampleType::hijing);
	// TEST sample: reads the overlay test sample AND uses the pp beam only (only beam produced).
	// ONE switch -- it also selects the input dir. See FullSimSampleType.h.
	py.isTestSample = true;
	py.overlay_pbpb_year = 23;   // label hijing_overlay_pbpb23, matching INDIR above
	// AMI PROVENANCE: the overlay is built on the pp-beam evgen DSIDs. AMI files are keyed by
	// beam+slice only, so declare them; InitInputFullsim throws if another production is read.
	py.expected_ami_dsids = {802776, 802777, 802778, 802779, 802780, 802781};
	// Single-slice DIAGNOSTIC run: these input dirs hold pTH8_14 ONLY. Every output is a
	// ratio, so a missing slice is intentional and harmless here.
	py.allow_missing_slices = true;
	py.fullsim_input_dir_override = "${INDIR}";
	py.fill_kn_trees_fullsim = true;
	py.output_single_muon_tree = ${SINGLE};
	py.extra_output_suffix = "${SUFFIX}";
	py.Run();

	.q;

EOF
