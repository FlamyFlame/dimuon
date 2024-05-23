# setupATLAS

ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh ""

lsetup "views LCG_102 x86_64-centos7-gcc11-opt"

export LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH
# export PATH=$PYTHIA8/bin:$PATH
