#!/bin/bash

# Run Pythia fullsim HIJING overlay: single-muon tree output for reco efficiency sanity check.
# Usage: ./run_pythia_fullsim_overlay_single_muon.sh

cd "$(dirname "$0")"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

: "${OVERLAY_YEAR:=24}"

root -b -l << EOF
	.L PythiaAnalysisClasses.h+

	PythiaFullSimOverlayAnalysis py(FullSimSampleType::hijing);
	// TEST sample: reads the overlay test sample AND uses the pp beam only (only beam produced).
	// ONE switch -- it also selects the input dir. See FullSimSampleType.h.
	py.isTestSample = true;
	// HIJING-overlay conditions year: 24 (default) = pythia_fullsim_hijing_overlay_test_sample/
	// (pbpb24, r17864), 23 = ..._pbpb23/ (r17618). Env OVERLAY_YEAR; drives dir AND label.
	py.overlay_pbpb_year = ${OVERLAY_YEAR};
	// AMI PROVENANCE: the overlay is built on the pp-beam evgen DSIDs. AMI files are keyed by
	// beam+slice only, so declare them; InitInputFullsim throws if another production is read.
	py.expected_ami_dsids = {802776, 802777, 802778, 802779, 802780, 802781};
	// SINGLE-SLICE pbpb24 TEST sample (r17864, pTH125_300 only; hijing_overlay_pbpb24_test_sample_skim.md):
	// a missing pT-hat slice is fatal by default because it biases the sigma-weighted combination.
	// TEMPORARY, diagnostic only: allow it for year 24 until the full 6-slice pbpb24 production
	// exists -- then DELETE this line so the guard is strict again. The 6-slice pbpb23 sample stays strict.
	py.allow_missing_slices = (${OVERLAY_YEAR} == 24);
	py.fill_kn_trees_fullsim = true;
	py.output_single_muon_tree = true;
	py.extra_output_suffix = "_single_muon";
	py.Run();

	.q;

EOF
