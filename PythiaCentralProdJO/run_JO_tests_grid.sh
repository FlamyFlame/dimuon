#!/usr/bin/env bash
# make_and_submit.sh ─ generate Condor wrappers + submit
# ──────────────────────────────────────────────────────
set -euo pipefail

# setupATLAS
# asetup AthGeneration,23.6.47 # set up athena package locally

joDir="JOs_5.36TeV"   # directory where the 24 JO files live

# ---- definitions -----------------------------------------------------------
# pthLabels=(5_8 8_14 14_24 24_40 40_70 70_125 125_300)
# jobCounts=(10 5 5 1 1 1 1) # jobs per slice
# jobCounts=(2 2 2 2 1 1 1) # jobs per slice
pthLabels=(5_8)
jobCounts=(10)
isospins=(pp pn np nn)

# ---- create + submit -------------------------------------------------------
for pIdx in "${!pthLabels[@]}"; do
  pTag=${pthLabels[$pIdx]}
  nJobs=${jobCounts[$pIdx]}

  for iIdx in "${!isospins[@]}"; do
    iso=${isospins[$iIdx]}
    cfg="${iso}_${pTag}"                 # e.g. pp_25_60
    workDir="JO_test_${cfg}" # name for both condor-script & log-file directories in /afs, and for JO & output-EVNT-file directories in /eos

    joFile="mc.Py8EG_A14_5360_${iso}_hQCD_DiMu_pTH${pTag}.py"

    for jobIdx in $(seq 1 "$nJobs"); do
      seed=$(( (${pIdx}+1)*10000 + ${iIdx}*1000 + jobIdx ))
      
      rm -rf ${workDir}/job${jobIdx}
      mkdir -p ${workDir}/job${jobIdx}

      cp ${joDir}/${joFile} ${workDir}/job${jobIdx}/
      cd ${workDir}/job${jobIdx}

      # ---------- executable ---------------------------------------------------
      cat > "./run_JO_test.sh" <<EOF
#!/usr/bin/env bash

start_time=\$(date +%s)

Gen_tf.py --randomSeed=${seed} \\
          --jobConfig=\$PWD \\
          --outputEVNTFile="Pythia.${seed}.EVNT.pool.root" > output.txt

end_time=\$(date +%s)
elapsed=\$((end_time - start_time))
echo "Job completed in \$elapsed seconds" >> output.txt

EOF
      chmod +x run_JO_test.sh

      # ---------- run grid job -----------------------------------------------------
      prun --exec "./run_JO_test.sh" \
           --useAthenaPackage \
           --outDS user.yuhang.PyJO.dimuon.May25.${iso}.pTHat${pTag}GeV.job${jobIdx}.v7 \
           --outputs Pythia.${seed}.EVNT.pool.root,output.txt,eventLoopHeartBeat.txt

      cd ../..
    done
  done
done
