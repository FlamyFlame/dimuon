#!/usr/bin/env bash
# make_and_submit.sh ─ generate Condor wrappers + submit
# ──────────────────────────────────────────────────────
set -euo pipefail

joDir="/eos/user/y/yuhang/dimuon/PythiaCentralProdJO/JOs"   # directory where the 24 JO files live

# ---- definitions -----------------------------------------------------------
pthLabels=(5_8 8_14 14_24 24_40 40_70 70_125 125_300)
jobCounts=(1 1 1 1 1 1 1) # jobs per slice
isospins=(pp pn np nn)

# ---- create + submit -------------------------------------------------------
for pIdx in "${!pthLabels[@]}"; do
  pTag=${pthLabels[$pIdx]}
  nJobs=${jobCounts[$pIdx]}

  for iIdx in "${!isospins[@]}"; do
    iso=${isospins[$iIdx]}
    cfg="${iso}_${pTag}"                 # e.g. pp_25_60
    workDir="JO_test_${cfg}" # name for both condor-script & log-file directories in /afs, and for JO & output-EVNT-file directories in /eos
    tmpDir="/tmp/yuhang/${workDir}"
    outDir="/eos/user/y/yuhang/python_jo_test_area/${workDir}"
    mkdir -p "${workDir}/log"
    mkdir -p "${outDir}"
    chmod a+w "${outDir}"
    joFile="mc.Py8EG_A14_${iso}_HardQCD_DiMu_pTH${pTag}.py"

    # ---------- executable ---------------------------------------------------
    cat > "${workDir}/run_JO_test_${cfg}.sh" <<EOF
#!/usr/bin/env bash
# \$1 = Condor process number (0– …)
jobIdx=\$(( \$1 + 1 ))                 # 1‑based inside slice
seed=\$(( (${pIdx}+1)*10000 + ${iIdx}*1000 + jobIdx ))

mkdir -p "${tmpDir}"
chmod a+w "${tmpDir}"
cp "${joDir}/${joFile}" "${tmpDir}/"
cd ${tmpDir}

mkdir -p job\${jobIdx}
chmod a+w job\${jobIdx}

cp "${joFile}" "job\${jobIdx}/"
cd job\${jobIdx}

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source \$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
asetup 23.6.47,AthGeneration

start_time=\$(date +%s)

Gen_tf.py --randomSeed=\${seed} \\
          --jobConfig=\$PWD \\
          --outputEVNTFile="Pythia.\${seed}.EVNT.pool.root"

end_time=\$(date +%s)
elapsed=\$((end_time - start_time))
echo "Job completed in \$elapsed seconds"

mv Pythia.\${seed}.EVNT.pool.root ${outDir}/
mv eventLoopHeartBeat.txt ${outDir}/
mv log.generate ${outDir}/

EOF
    chmod +x "${workDir}/run_JO_test_${cfg}.sh"

    # ---------- condor submit -----------------------------------------------
    cat > "${workDir}/run_JO_test_${cfg}.sub" <<EOF
universe        = vanilla
executable      = run_JO_test_${cfg}.sh
arguments       = \$(Process)
output          = log/out_\$(Process).txt
error           = log/err_\$(Process).txt
log             = log/condor_\$(Process).log
getenv          = True
request_memory  = 4GB
+JobFlavour     = "tomorrow"
queue ${nJobs}
EOF

    # ---------- queue it -----------------------------------------------------
    ( cd "${workDir}" && condor_submit "run_JO_test_${cfg}.sub" )
  done
done
