#!/bin/bash

export TargetDir="$PWD"/condorRun

rm inputdata.txt

if [ -d ${TargetDir} ]; then
  rm -rf ${TargetDir}/OutDir*
else
  mkdir ${TargetDir}
fi


#for file in /atlasgpfs01/usatlas/data/bseidlit/trees/pp_5TeV_w_tracksClusters_HighMu/user.b*.root/*
#for file in /atlasgpfs01/usatlas/scratch/bseidlit/pp_w_tracksAndClus/user.bseidlit.ppDecorrelation_14.0034*.r11215_p3764_ANALYSIS.root/*
#for file in /usatlas/scratch/bseidlit/pp_13TeV_PYTHIA/user.bseidlit.ppDecorrelation_mb_pythia_01.361203.e3639_s2601_s2132_r6616_r6270_tid05411175_00_ANALYSIS.root/* 
#for file in /sphenix/u/bseidlitz/work/trees/pp_5TeV_tracksClusters/*/*
for file in /atlasgpfs01/usatlas/data/bseidlit/trees/pp_13TeV_w_trackClusters_moreClus/*/*
do
cat >>inputdata.txt<< EOF
$file
EOF
done


j=1000
tot_files=$( cat inputdata.txt | wc -l )
echo "total files: $tot_files"
rem=$(( $tot_files%$j ))
files_per_job=$(( $tot_files/$j ))
njob=$j
if [ $rem -ne 0 ]; then
  files_per_job=$(( $files_per_job+1 ))
fi
rem2=$(( $tot_files%$files_per_job ))
njob=$(( $tot_files/$files_per_job ))
if [ $rem2 -ne 0 ]; then
  njob=$(( ($tot_files/$files_per_job)+1 ))
fi
echo "files per job: $files_per_job"
echo "njob: $njob"

for((i=0;i<$njob;i++));
do

  mkdir ${TargetDir}/OutDir$i
  export WorkDir="${TargetDir}/OutDir$i"
  echo "WorkDir:" ${WorkDir}
  start_file=$(( $i*$files_per_job+1 ))
  end_file=$(( $start_file+$files_per_job-1 ))
  echo "start file: $start_file   end file: $end_file"

  sed -n $start_file\,${end_file}p inputdata.txt > tmp.txt
  mv tmp.txt ${WorkDir}/inputdata.txt

  pushd ${WorkDir}

  cp "$PWD"/../../CondorRun.sh CondorRunTC$i.sh
  cp -v "$PWD"/../../runMaker .
  cp "$PWD"/../../trkEff_nominal.root .

  cat >>ff.sub<< EOF
+JobFlavour                   = "workday"
transfer_input_files          = ${WorkDir}/runMaker,${WorkDir}/inputdata.txt,${WorkDir}/CondorRunTC$i.sh,${WorkDir}/trkEff_nominal.root
Executable                    = CondorRunTC$i.sh
Universe                      = vanilla
Notification                  = Never
GetEnv                        = True
Priority                      = +20
Output                        = test.out
Error                         = test.err
Log                           = test.log
Notify_user                   = blair.daniel.seidlitz@cern.ch
accounting_group              = group_atlas.boulder

Queue
EOF

  condor_submit ff.sub
  popd
done
