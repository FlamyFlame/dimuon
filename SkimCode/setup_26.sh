export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh

# AthAnalysis 25.2.90 -- the FIRST release whose TrigConfData knows the gFEX
# EnergyThreshold flavour "gRISTRETTO" that appears in the 2026 L1 menu.  25.2.89
# (used for the 2023/24/25 skims) aborts on the first 2026 event with
#   "Flavour gRISTRETTO for EnergyThreshold algorithm not recongnised!"
# 25.2.90 keeps the same gcc14 / LCG_108a_ATLAS_9 platform as 25.2.89, so it is the
# smallest possible step away from the release the other years were skimmed with.
#
# Same cmake pin as setup_25.sh: ALRB's current cmake 4.x is incompatible with this
# LCG's FindPython module.
_CMAKE329=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Cmake/3.29.5/Linux-x86_64/bin
export PATH=$_CMAKE329:$PATH

mkdir -p build_26
cd build_26
acmSetup --sourcedir=../source AthAnalysis,25.2.90
export PATH=$_CMAKE329:$PATH

acm compile
cd ..
