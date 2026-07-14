#!/bin/bash

# Pythia pp24-conditions fullsim -- the FULL sample ("_pdf" production, DSIDs 803015-803020).
#
# Reads the LOCALGROUPDISK symlink farm (the ~280 GB of NTUP is never hadded onto GPFS:
# see Analysis/docs/tracking/pythia_fullsim_pp24_full_sample_skim.md, Design Decision D1).
# Each pT-hat slice is N unmerged grid outputs exposed as "...NTUP.partNN.root"; the reader
# TChain-globs them.
#
# THREE settings distinguish this from the TEST-sample runs -- all three are mandatory:
#   1. fullsim_input_dir_override -> the farm (not pythia_fullsim_test_sample/)
#   2. ami_info_dir_override + expected_ami_dsids -> the "_pdf" production's OWN cross-sections.
#      AMI files are keyed by beam+slice ONLY, so without this the sample would silently be
#      weighted with the TEST sample's sigma*eff (pTH8_14: 23.78 nb vs the correct 36.74 nb --
#      a 55% error that VARIES BY SLICE and therefore does NOT cancel in ratios).
#      expected_ami_dsids makes InitInputFullsim THROW if the wrong production is picked up.
#   3. extra_output_suffix -> so this never clobbers the test sample's
#      muon_pairs_pythia_fullsim_pp24*.root
#
# Isospin: pp CONDITIONS simulate pp collisions -> the pp beam alone, isospin weight 1.
# That is the DEFAULT for FullSimSampleType::pp, so this script sets nothing (unlike the TEST
# sample, which was produced with 4 beams by mistake and must opt in with isTestSample=true).
#
# Usage: ./run_pythia_fullsim_full_sample.sh
# Prereq: the farm exists -- SkimCode/scripts/fullsim_pp24_full_to_lgd.sh

cd "$(dirname "$0")"

FARM_DIR="/usatlas/u/yuhanguo/usatlasdata/pythia_fullsim_full_sample/"
AMI_DIR="${FARM_DIR}ami_info/"

if ! ls "${FARM_DIR}"Pythia_5p36TeV_pp_hQCD_DiMu_pTH*.FullSimPP24.NTUP.part*.root >/dev/null 2>&1; then
    echo "ERROR: no symlink farm in ${FARM_DIR}"
    echo "       Run SkimCode/scripts/fullsim_pp24_full_to_lgd.sh first."
    exit 1
fi

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

root -b -l << EOF
	.L PythiaAnalysisClasses.h+

	PythiaFullSimAnalysis py(0);
	// isTestSample defaults to FALSE = the FULL production => input dir, AMI dir and the
	// pp-beam-only/weight-1 isospin all follow from that one switch. Do NOT set it true.
	py.fullsim_input_dir_override = "${FARM_DIR}";
	py.ami_info_dir_override      = "${AMI_DIR}";
	py.expected_ami_dsids         = {803015, 803016, 803017, 803018, 803019, 803020};
	py.extra_output_suffix        = "_full";
	py.fill_kn_trees_fullsim      = true;
	py.Run();
	.q;
EOF
