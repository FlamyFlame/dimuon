#!/bin/bash
# Diagnostic: print TTree::Print() for ZDC/FCal branches and scan values.
# Run interactively: bash run_pbpb_24_diagnose.sh

cd /gpfs/mnt/atlasgpfs01/usatlas/workarea/yuhanguo/dimuon_codes/Analysis/NTupleProcessingCode

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

lsetup "views LCG_107a_ATLAS_2 x86_64-el9-gcc13-opt"

root -b -q test_check_skim_branches.cxx
