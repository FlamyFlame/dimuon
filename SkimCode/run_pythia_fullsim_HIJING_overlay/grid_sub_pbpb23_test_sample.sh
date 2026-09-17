#!/usr/bin/env bash

# Pythia8 HardQCD pp -> DiMu Fullsim, reweighted HIJING overlay (2024 HI conditions)
# mc23_5p36TeV AOD containers: r17618 reco over r15970 HIJING merge.
# 6 datasets, 10 files each (~83-85 GB / 10k events per dataset).
#
# r17618 steering is "doRDO_TRIG" "doTRIGtoALL" with Trigger.AODEDMSet='AODFULL' and
# HLT menu Dev_HI_run3_v1 => trigger decisions + HLT muon containers are present, so
# this skim runs with the trigger ENABLED (July2026 and later).
#
# Output file name is driven by TRIGRATES_OUTPUT_FILE and matches --extOutFile:
#   Pythia_5p36TeV_<beam>_hQCD_DiMu_pTH<lo>_<hi>.<OUT_TAG>.NTUP.root
#
# NTUP NAMING (user decision 2026-09-17): OUT_TAG = FullSimHIJINGOverlayPbPb<yy>.<cfg>,
# <cfg> = vtxz<z>mm_b<lo>_<hi>fm = the production configuration (pinned vertex z in mm,
# HIJING impact-parameter interval from the HITS "ip" token) -- full description in
# grid_sub_pbpb24_test_sample.sh. This sample is ONE point: r17618 vtx (-0.6,-0.4,-3.3) mm,
# ip0_5 -> FullSimHIJINGOverlayPbPb23.vtxz-3_3mm_b0_5fm.
# History: the July2026.v2 tasks were submitted with the legacy tag FullSimHIJINGOverlayPP24
# ("PP24" was wrong -- Pb+Pb 2023 conditions); the downloaded NTUPs were mv'd to the new
# name on 2026-09-17. The GRID datasets (and the .bak_20260709 LGD copies) keep the old name.
#
# Run from: SkimCode/run_pythia_fullsim_HIJING_overlay
# Setup:    source ../setup_25.sh && lsetup panda

_cfg_tag() { echo "vtxz${1//./_}mm_b${2}fm"; }   # -3.3 0_5 -> vtxz-3_3mm_b0_5fm

COND_YEAR=23
# VER_TAG = the version of the task that produced the files on disk (history record; this
# sample is not resubmitted). A resubmission would bump it, as grid_sub_pbpb24_test_sample.sh does.
VER_TAG=July2026.v2

_submit() {
  local sample="$1"     # e.g. pp_hQCD_DiMu_pTH125_300
  local inDS="$2"       # rucio scope:container
  local vtxz="$3"       # pinned vertex z [mm]
  local ip="$4"         # HIJING HITS impact-parameter token
  local out_tag="FullSimHIJINGOverlayPbPb${COND_YEAR}.$(_cfg_tag "$vtxz" "$ip")"
  local fname="Pythia_5p36TeV_${sample}.${out_tag}.NTUP.root"
  local outDS="user.yuhang.NTUP.Pythia_5p36TeV_${sample}.${out_tag}.${VER_TAG}."
  pathena --trf "TRIGRATES_RUNMODE=ppmcfullsim_hioverlay24 TRIGRATES_OUTPUT_FILE=${fname} athena.py TrigRates_CA.py --filesInput=%IN --evtMax=%MAXEVENTS" \
          --inDS "${inDS}" \
          --outDS "${outDS}" \
          --mergeOutput \
          --nFilesPerJob 2 \
          --extOutFile "${fname}" \
          --excludeFile=clean.sh \
          --excludeFile=vim_backup
}

_submit pp_hQCD_DiMu_pTH125_300 mc23_5p36TeV:mc23_5p36TeV.802776.Py8EG_A14_pp_hQCD_DiMu_pTH125_300.merge.AOD.e8599_s4614_r17618_r15970 -3.3 0_5
_submit pp_hQCD_DiMu_pTH14_24   mc23_5p36TeV:mc23_5p36TeV.802777.Py8EG_A14_pp_hQCD_DiMu_pTH14_24.merge.AOD.e8599_s4614_r17618_r15970 -3.3 0_5
_submit pp_hQCD_DiMu_pTH24_40   mc23_5p36TeV:mc23_5p36TeV.802778.Py8EG_A14_pp_hQCD_DiMu_pTH24_40.merge.AOD.e8599_s4614_r17618_r15970 -3.3 0_5
_submit pp_hQCD_DiMu_pTH40_70   mc23_5p36TeV:mc23_5p36TeV.802779.Py8EG_A14_pp_hQCD_DiMu_pTH40_70.merge.AOD.e8599_s4614_r17618_r15970 -3.3 0_5
_submit pp_hQCD_DiMu_pTH70_125  mc23_5p36TeV:mc23_5p36TeV.802780.Py8EG_A14_pp_hQCD_DiMu_pTH70_125.merge.AOD.e8599_s4614_r17618_r15970 -3.3 0_5
_submit pp_hQCD_DiMu_pTH8_14    mc23_5p36TeV:mc23_5p36TeV.802781.Py8EG_A14_pp_hQCD_DiMu_pTH8_14.merge.AOD.e8599_s4614_r17618_r15970 -3.3 0_5
