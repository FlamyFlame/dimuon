#!/bin/bash

# Pythia pp24 fullsim FULL sample -- single-muon tree (denominator for the single-muon reco efficiency).
#
# FULL SAMPLE (the "_pdf" production, DSIDs 803015-803020) -- reads the LOCALGROUPDISK symlink
# farm; the ~280 GB of NTUP is never hadded onto GPFS (tracking doc
# pythia_fullsim_pp24_full_sample_skim.md, Design Decision D1).
#
# isTestSample is left FALSE (the default) => input dir, AMI dir (the proton-PDF evgen's
# ami_info_PDF/, FullSimSampleType.h) and the pp-beam-only /
# isospin-weight-1 treatment ALL follow from that one switch (FullSimSampleType.h).
# expected_ami_dsids is the BLOCKING provenance guard: AMI files are keyed by beam+slice only,
# so without it this run could be silently weighted with the TEST sample's cross-sections
# (slice-dependent => cancels NOWHERE). See Analysis/docs/ami_weights.md.
#
# SUFFIX RULE: "_full" is always the LAST suffix, so every full-sample product ends in _full and
# consumers can simply append the sample suffix.
# Output: muon_pairs_pythia_fullsim_pp24_no_data_resonance_cuts_single_muon_full.root
#
# Usage: ./run_pythia_fullsim_single_muon_full_sample.sh

cd "$(dirname "$0")"

FARM_DIR="/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/"
# Smoke test: NEVENTS_MAX=2000 ./this_script.sh   (-1 = all events, the default)
NEVENTS_MAX="${NEVENTS_MAX:--1}"

if ! ls "${FARM_DIR}"Pythia_5p36TeV_pp_hQCD_DiMu_pTH*.FullSimPP24.NTUP.part*.root >/dev/null 2>&1; then
    echo "ERROR: no symlink farm in ${FARM_DIR}"
    echo "       Run SkimCode/scripts/fullsim_pp24_full_to_lgd.sh first."
    exit 1
fi

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

root -b -l << ROOTEOF
	.L PythiaAnalysisClasses.h+

	PythiaFullSimAnalysis py(0);
	// isTestSample stays FALSE (default) = the FULL production. Do NOT set it true.
	py.fullsim_input_dir_override = "${FARM_DIR}";
	py.expected_ami_dsids         = {803015, 803016, 803017, 803018, 803019, 803020};
	py.fill_kn_trees_fullsim      = true;
	py.output_single_muon_tree    = true;
	py.extra_output_suffix        = "_single_muon_full";
	py.nevents_max                = ${NEVENTS_MAX};
	py.Run();
	.q;
ROOTEOF
