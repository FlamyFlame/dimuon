export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

cd build_25
acmSetup --sourcedir=../source_25 AthAnalysis,25.2.35
cd ..
