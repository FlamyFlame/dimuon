export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

# AthAnalysis 25.2.89 uses LCG_108a_ATLAS_9 (gcc14).  ALRB's "current" cmake
# (4.2.1) is incompatible with that LCG's FindPython module: cmake reports
# "Cannot run the interpreter python3.12" and aborts configuration.
# Override to cmake 3.29.5 BEFORE acmSetup (so fresh builds configure with it)
# and AGAIN after (in case acmSetup re-prepends 4.x).
_CMAKE329=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Cmake/3.29.5/Linux-x86_64/bin
export PATH=$_CMAKE329:$PATH

cd build_25
acmSetup --sourcedir=../source AthAnalysis,25.2.89
export PATH=$_CMAKE329:$PATH  # re-assert in case acmSetup re-added 4.x

# Run acm compile so the local InstallArea (preinstall) is always up-to-date.
# pathena's cpack calls `cmake --build . --target preinstall` internally;
# if that target is stale the cpack subprocess fails.
acm compile
cd ..
