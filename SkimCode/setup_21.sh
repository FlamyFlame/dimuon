export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

cd build_21
acmSetup --sourcedir=../source_21 AthAnalysis,21.2.200
cd ..
