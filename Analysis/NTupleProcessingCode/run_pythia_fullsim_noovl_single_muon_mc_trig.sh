#!/bin/bash

# r17663 NO-OVERLAY diagnostic sample, trigger propagation (store_mc_trigger): SINGLE-MUON tree.
#
# Purpose (mc_trigger_efficiency.md R8 / round-4 Autonomy Contract): r17663 = r17618 minus
# HIJING (PileUp=False), everything else identical (Athena 24.0.58, OFLCOND-MC23-SDR-RUN3-05,
# ConditionsRunNumber=460000, L1 MC_HI_run3_v1). Single-variable discriminator for the
# forward-endcap L1 anomaly: does a QUIET event with HI-conditions reco show the pp24-r16578
# saturation (~0.9) or the r17618 turn-on (~0.5) in q.eta in (-2.4,-2), pT 4-6?
#
# Sample identity is entirely separate from pp24 and the overlay: own input dir, own file
# tag (FullSimHIJINGOverlayPP24_r17663), own output label (r17663_no_overlay) -> NOTHING
# existing is overwritten.
#
# NOT an overlay: no HIJING, no centrality -> the pp analysis class (no overlay Extras).
# Usage: ./run_pythia_fullsim_noovl_single_muon_mc_trig.sh

cd "$(dirname "$0")"

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh --quiet
lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

root -b -l << EOF
	.L PythiaAnalysisClasses.h+

	// sample_type = noovl -> input dir / file tag / label all switch together.
	PythiaFullSimAnalysis py(0, false, false, FullSimSampleType::noovl);
	// Only ONE pT-hat slice (pTH8_14, DSID 802781) was produced. Every other slice is
	// legitimately absent, so the missing-slice guard must be lifted. This is a
	// DIAGNOSTIC sample: with a single slice the sigma*eff weight is one global constant
	// that cancels in every efficiency ratio -- no cross-section-weighted result is taken
	// from this run.
	py.allow_missing_slices = true;
	// AMI PROVENANCE: r17663 uses the SAME evgen dataset as the overlay/test pTH8_14 pp
	// slice (mc23_5p36TeV.802781 ... e8599; verified in AMI 2026-07-16: crossSection
	// 4.816e+06 nb, genFiltEff 4.938681e-06 -- identical, because the AMI weight lives on
	// the EVGEN dataset and the r-tag does not change it). Declare the DSID; the code
	// throws if the AMI file it reads carries any other datasetNumber.
	py.expected_ami_dsids = {802781};
	// isTestSample=true -> nPDF evgen -> AMI dir = the evgen's ami_info_nPDF/ (where the 802781
	// file lives; FullSimSampleType.h). The input dir is the same either way for this sample type.
	py.isTestSample = true;
	py.output_single_muon_tree = true;
	py.extra_output_suffix = "_single_muon";
	py.store_mc_trigger = true;
	py.Run();

	.q;

EOF
