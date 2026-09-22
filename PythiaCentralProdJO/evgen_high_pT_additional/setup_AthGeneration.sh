#!/bin/bash
# Setup ATLAS environment
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

cd build
acmSetup --sourcedir=../source/athena AthGeneration,23.6.68

# Run acm compile so the local InstallArea (preinstall) is always up-to-date.
# pathena's cpack calls `cmake --build . --target preinstall` internally;
# if that target is stale the cpack subprocess fails.
acm compile
cd ..