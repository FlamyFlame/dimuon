#!/usr/bin/env bash

# Pythia8 HardQCD pp -> DiMu Fullsim, HIJING overlay with **Pb+Pb 2024 conditions**
# (r17864 reco over r17855 HIJING merge; evgen e8613_e8586; sim s4684).
#   r17864: "2024 run conditions vtx(-.7,-.6,-71.2) for ATLHI-576, im_par=0-5
#            - bunch structure added to r17835"; preInclude Campaigns.MC23HeavyIons2024NoPileUp,
#            HIJING HITS mc23_5p36TeV.860250.Hijing_PbPb_UCC_Flow_JJFV6_ip0_5 (b = 0-5 fm),
#            steering doRDO_TRIG + doTRIGtoALL (trigger simulated, 2024 HI menu).
# TEST production: ONE DSID only (pTH125_300), 10 files / 10k events / 54.9 GB.
# The previous test sample (Pb+Pb 2023 conditions, r17618, 6 slices) is
# grid_sub_pbpb23_test_sample.sh; its NTUPs now live in
# ~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample_pbpb23/.
#
# AMI weights are the PYTHIA EVGEN's (docs/ami_weights.md, user ruling 2026-09-17): this
# sample is 802776 e8599 (nPDF evgen, registry table A); e8613 in its tag chain is the HIJING
# evgen, not a Pythia one. The AOD-level AMI record fetched 2026-09-15 (sigma 89.541 nb,
# genFiltEff 2.342314e-3) is kept as ~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/
# ami_info_AOD_record_20260915/ and is read by nothing.
#
# Output file name is driven by TRIGRATES_OUTPUT_FILE and matches --extOutFile:
#   Pythia_5p36TeV_<beam>_hQCD_DiMu_pTH<lo>_<hi>.<OUT_TAG>.NTUP.root
#
# NTUP NAMING (user decision 2026-09-17). OUT_TAG = FullSimHIJINGOverlayPbPb<yy>.<cfg> where
# <cfg> = vtxz<z>mm_b<lo>_<hi>fm names the production CONFIGURATION on top of the isospin
# beam and pT-hat slice already in the sample name: the pinned vertex z in mm ('.' -> '_',
# sign kept) and the HIJING impact-parameter interval (the HITS "ip" range, ip0_5 -> b0_5fm).
# Pb+Pb 2024 requests FIVE z points -71.2 / -28.5 / -4.8 / 18.9 / 61.6 mm (x = -0.7,
# y = -0.6 mm always) x FOUR b intervals 0-5 / 5-9 / 9-12 / 12+ fm (the open one takes its
# HITS token when it exists); each point is its own dataset and its own NTUP -> one
# _submit line per (beam, slice, z, b). The reader side is
# Analysis/MuonObjectsParamsAndHelpers/FullSimSampleType.h (FullSimSampleFileTag /
# FullSimOverlayConfigTag) -- keep the two in step.
# History: jediTaskID 52568863 (VER_TAG Sep2026.v1) was submitted with the legacy tag
# FullSimHIJINGOverlayPP24 ("PP24" was wrong -- the sample is Pb+Pb 2024); its downloaded
# NTUP was mv'd to the new name on 2026-09-17. The GRID dataset keeps the old name.
#
# Run from: SkimCode/run_pythia_fullsim_HIJING_overlay
# Setup:    source ../setup_25.sh && lsetup panda
# Monitor:  ../scripts/grid_monitor.sh --mode overlay <taskid>
#           (writes flat into ~/usatlasdata/pythia_fullsim_hijing_overlay_test_sample/)

# _cfg_tag <vtx z in mm, e.g. -71.2> <HITS ip range, e.g. 0_5>  ->  vtxz-71_2mm_b0_5fm
_cfg_tag() { echo "vtxz${1//./_}mm_b${2}fm"; }

COND_YEAR=24
# VER_TAG = skim-code / submission version, bumped on EVERY resubmission of the same input (a
# new OUT_TAG alone would already give a distinct outDS, but the version keeps the history
# linear: v1 = the legacy-tag task above, v2 = first submission under the new naming).
VER_TAG=Sep2026.v2

_submit() {
  local sample="$1"     # e.g. pp_hQCD_DiMu_pTH125_300
  local inDS="$2"       # rucio scope:container
  local vtxz="$3"       # pinned vertex z [mm], as in the r-tag description
  local ip="$4"         # HIJING HITS impact-parameter token, e.g. 0_5
  local out_tag="FullSimHIJINGOverlayPbPb${COND_YEAR}.$(_cfg_tag "$vtxz" "$ip")"
  local fname="Pythia_5p36TeV_${sample}.${out_tag}.NTUP.root"
  local outDS="user.yuhang.NTUP.Pythia_5p36TeV_${sample}.${out_tag}.${VER_TAG}."
  pathena --trf "TRIGRATES_RUNMODE=ppmcfullsim_hioverlay24 TRIGRATES_OUTPUT_FILE=${fname} athena.py TrigRates_CA.py --filesInput=%IN --evtMax=%MAXEVENTS" \
          --inDS "${inDS}" \
          --outDS "${outDS}" \
          --mergeOutput \
          --nFilesPerJob 2 \
          --extOutFile "${fname}" \
          --excludeFile="clean.sh,vim_backup,test_r17618,test_r17662,test_r17663,test_r17864,mc23_5p36TeV,\*.root,\*.log"
}

# r17864: vtx z = -71.2 mm, HITS ip0_5
_submit pp_hQCD_DiMu_pTH125_300 mc23_5p36TeV:mc23_5p36TeV.802776.Py8EG_A14_pp_hQCD_DiMu_pTH125_300.merge.AOD.e8613_e8586_s4684_r17864_r17855_tid52480369_00 -71.2 0_5
