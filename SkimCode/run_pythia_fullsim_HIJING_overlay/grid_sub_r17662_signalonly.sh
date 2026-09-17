#!/usr/bin/env bash
#
# Signal-only-truth (r17662) HIJING-overlay skim, pTH8_14, FULL statistics.
#
# Purpose: r17662 keeps ONLY the Pythia signal truth in the TruthParticleContainer
# (no HIJING truth), so the Pythia/HIJING barcode collision present in r17618 does
# NOT occur. This lets prob>0.5 matching work unambiguously — reproducing the Run 2
# note method (STARLIGHT+HIJING overlay, prob>0.5) which gives physical efficiencies.
#
# Input AOD: mc23_5p36TeV.802781 ... r17662_r15970 (10 files, 10000 events; verified
# via `rucio list-files`). The earlier 100-event signal-only skim came from
# run_test.sh's hardcoded --evtMax=100; here we process all events (--evtMax=%MAXEVENTS).
#
# Run from: SkimCode/run_pythia_fullsim_HIJING_overlay
# Setup:    source ../setup_25.sh && lsetup panda

# NTUP NAMING (2026-09-17): FullSimHIJINGOverlayPbPb<yy>_r17662.<cfg>, <cfg> = the production
# configuration vtxz<z>mm_b<lo>_<hi>fm (see grid_sub_pbpb24_test_sample.sh). r17662 = same
# conditions as r17618: Pb+Pb 2023, vtx (-0.6,-0.4,-3.3) mm, HITS ip0_5. The July2026.v1 grid
# dataset was made with the legacy tag FullSimHIJINGOverlayPP24_r17662; the downloaded NTUP
# was mv'd to the new name on 2026-09-17.
OUT_TAG=FullSimHIJINGOverlayPbPb23_r17662.vtxz-3_3mm_b0_5fm
VER_TAG=July2026.v1

sample=pp_hQCD_DiMu_pTH8_14
inDS=mc23_5p36TeV:mc23_5p36TeV.802781.Py8EG_A14_pp_hQCD_DiMu_pTH8_14.merge.AOD.e8599_s4614_r17662_r15970_tid50427580_00
fname="Pythia_5p36TeV_${sample}.${OUT_TAG}.NTUP.root"
outDS="user.yuhang.NTUP.Pythia_5p36TeV_${sample}.${OUT_TAG}.${VER_TAG}."

pathena --trf "TRIGRATES_RUNMODE=ppmcfullsim_hioverlay24 TRIGRATES_OUTPUT_FILE=${fname} athena.py TrigRates_CA.py --filesInput=%IN --evtMax=%MAXEVENTS" \
        --inDS "${inDS}" \
        --outDS "${outDS}" \
        --mergeOutput \
        --nFilesPerJob 2 \
        --extOutFile "${fname}" \
        --excludeFile=clean.sh \
        --excludeFile=vim_backup
