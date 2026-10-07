#!/bin/bash
# Local EVNT -> TRUTH0 -> SkimCode NTUP test, mirroring grid/chain.py of evgen_high_pT_additional.
# Usage: run_chain_test.sh <EVNT file> <outdir>
EVNT=$1; OUT=$2; mkdir -p $OUT/deriv $OUT/skim
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
( cd $OUT/deriv && set -- && source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh > /dev/null && asetup Athena,25.0.72 > /dev/null && \
  Derivation_tf.py --inputEVNTFile $EVNT --outputDAODFile output.pool.root --formats TRUTH0 > deriv_tf.out 2>&1; echo "DERIV_EXIT $?" >> deriv_tf.out )
SK=/usatlas/u/yuhanguo/workarea/dimuon_codes/SkimCode
( set -- && cd $SK && source $SK/setup_25.sh > /dev/null 2>&1 && cd $OUT/skim && ln -sf ../deriv/DAOD_TRUTH0.output.pool.root DAOD_TRUTH0.DAOD_TRUTH0.pool.root && \
  athena $SK/run_evgen_high_pT_additional/TruthTrigRates.py > skim.log 2>&1; echo "SKIM_EXIT $?" >> skim.log )
echo CHAIN_DONE
