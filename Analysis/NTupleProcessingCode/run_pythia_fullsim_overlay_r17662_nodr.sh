#!/bin/bash
# r17662 (signal-only-truth) overlay NTP, pTH8_14, FULL 10k events.
# Runs the DEFAULT procedure: pure prob>0.5 barcode matching, NO ad-hoc dR
# fallback (use_dr_fallback defaults false).  To enable the dR fallback (only
# meaningful for the barcode-collision r17618 sample), set py.use_dr_fallback=true.
#
# NOTE: do NOT try to control the dR fallback via pythia_only_barcode_cache —
# InitParamsExtra() re-sets it to true inside Run(), overriding any script value;
# and it only gates the ancestor-tracing cache, not the dR fallback.
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
	// DEFAULT: use_dr_fallback=false (no dR), use_geometric_matching not set -> pure prob>0.5
	py.Run();

	.q;

EOF
