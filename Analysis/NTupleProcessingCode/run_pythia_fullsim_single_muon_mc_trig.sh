#!/bin/bash

# Run Pythia pp24 fullsim with trigger propagation (store_mc_trigger):
# single-muon tree for the MC-based single-muon mu4 efficiency (Step 1,
# mc_trigger_efficiency.md §3.1). Gate = reco-matched, loose reco fiducial
# (exact pT>4, |eta|<2.4 + WP applied in RDF). Output suffix: _mc_trig_single_muon.
# Usage: ./run_pythia_fullsim_single_muon_mc_trig.sh

cd "$(dirname "$0")"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

root -b -l << EOF
	.L PythiaAnalysisClasses.h+

	PythiaFullSimAnalysis py(0);
	// TEST sample: reads pythia_fullsim_test_sample/ AND uses the 4-isospin-beam average
	// (that sample was produced with 4 beams by mistake). ONE switch -- see FullSimSampleType.h.
	py.isTestSample = true;
	// AMI PROVENANCE: AMI files are keyed by beam+slice only, so nothing otherwise stops
	// this run from being weighted with ANOTHER production's cross-sections. Declare the
	// TEST sample's DSIDs; InitInputFullsim throws on a mismatch.
	py.expected_ami_dsids = {802758, 802759, 802760, 802761, 802762, 802763,
                                802764, 802765, 802766, 802767, 802768, 802769,
                                802770, 802771, 802772, 802773, 802774, 802775,
                                802776, 802777, 802778, 802779, 802780, 802781};
	py.fill_kn_trees_fullsim = true;
	py.store_mc_trigger = true;
	py.output_single_muon_tree = true;
	py.extra_output_suffix = "_single_muon";
	py.Run();

	.q;

EOF
