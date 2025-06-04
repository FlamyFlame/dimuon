#!/bin/bash
# program idea: runs pythia simulation jobs with the job-name as input argument
# specifies code + output-file directories and job types (the nucleon-nucleon combinations: pp, pn, np, nn)
# compiles the latest .cc code, makes log + output-file directories for each job
# calculate #jobs for each (job type, kinematic range) combination 
# makes (kinematic range, job type)-specific logfile + outputfile sub-directories in the output area
# then for each job, copies over the executable, writes condor-submission script + bash script (sets up environment + runs the executable) , and condor_submit the condor script


# program requires a non-empty string as its first (and only) argument
if [ ${#1} -eq 0 ]
then
  echo "Must provide a non-empty string as argument!"
  echo "Suggestion: DATE_hardqcd_run# to avoid overwriting previous results"
  echo "e.g, 0315_hardqcd_run2"
  exit
fi

job_name=${1}
job_types=("pp" "pn" "np" "nn") # important: NO SPACE around "="
# njobs=(40 60 60 90) # array of number of jobs

PPJOBS=4
PNJOBS=6
NPJOBS=6
NNJOBS=9

subm_script_dir=$(pwd) # current directory where the submission scripts are

codedir=/usatlas/u/yuhanguo/workarea/dimuon_codes/PythiaSimulationCode
code=RunHardQCD_single_muon # The executable you want to run copies of

outputpath=/usatlas/u/yuhanguo/usatlasdata/pythia/${job_name} # Where all the results will end up
output=pytreeEff.root          # The output root tree name in each folder. Temporary to be moved to outputpath

rm -f $codedir/$code # Get rid of an old executable, if it exists
# g++ $codedir/$code.cc -o $codedir/$code -pedantic -W -Wall -Wshadow -fPIC -I$PYTHIA8_DIR/include -L$PYTHIA8_DIR/lib -lpythia8 -ldl -lRecursiveTools `root-config --libs --cflags` `/usr/local/bin/fastjet-config --cxxflags --libs --plugins` # Compile code with all necessary options. Current options are Pythia-8, fastjet, and soft-drop
g++ $codedir/$code.cc -o $codedir/$code -pedantic -W -Wall -Wshadow -fPIC -I$PYTHIA8/include -L$PYTHIA8/lib -lpythia8 -ldl `root-config --libs --cflags` `fastjet-config --cxxflags --libs --plugins` # Compile code with all 

rm -rf $outputpath    # If the output folder already existed, trash it.
mkdir -p $outputpath  # Make the output folder

mkdir -p $outputpath/jobs_${job_name}
cd $outputpath/jobs_${job_name}

echo "submitting job titled "${job_name}
file=c_job${job_name}.sub        # Definition of the file to contain the condor job instructions

rm -f $file           # If the condor file already existed, delete it.

echo "Universe        = vanilla" >>$file  # Tells condor the type of jobs it is queueing. Vanilla is standard
echo "GetEnv          = true" >>$file     # Tells condor to carry the current executing environment to new node
echo "Executable      = \$(Initialdir)/run.sh" >>$file    # The file condor will run as the job
echo "Output          = \$(Initialdir)/log/c_job.sub.out" >>$file  # Output location for the job's plaintext IO
echo "Error           = \$(Initialdir)/log/c_job.sub.err" >>$file  # Output location for the job's error IO 
echo "Log             = \$(Initialdir)/log/c_job.sub.log" >>$file  # Output location for condor's exec info

# Now, we will loop over all the jobs to set up the folder for its execution and add it to the condor job file
for k in $(seq 0 6); do # kinematic range

  mkdir -p $outputpath/k${k}

  if [ $k -eq 0 ]; then factor=4
  # elif [[ "$k" -eq 1 || "$k" -eq 2 ]]; then factor=40
  # elif [ $k -eq 3 ]; then factor=4
  # else factor=1
  else factor=0
  fi

  njobs=($[ PPJOBS * factor ] $[ PNJOBS * factor ] $[ NPJOBS * factor ] $[ NNJOBS * factor ])

  for itype in ${!njobs[@]}; do # loop over the indices of the njobs (and job_types) array
    njob=${njobs[itype]} # the itype-th element in the array njobs - number of jobs to generate
    job_type=${job_types[itype]} # the itype-th element in the array job_types - the current job type
  

    mkdir -p $outputpath/k${k}/${job_type}

    mkdir -p ${job_type}_k${k}_jobs_${job_name}
    # if [ ${k} -eq 0 ]; then
    #   echo "generating "${njob}" jobs in each kinematic range for job type "${job_type} # for sanity check
    # fi

    for job in $(seq 1 $njob); do # the itype-th job for a given beam setup (job type) and given kinematic range
      rm -rf ${job_type}_k${k}_jobs_${job_name}/job$job         # If the job folder already existed, delete it
      mkdir ${job_type}_k${k}_jobs_${job_name}/job$job          # Make a new job folder
      cp $codedir/$code ${job_type}_k${k}_jobs_${job_name}/job$job  # Copy the code executable into that folder
      mkdir ${job_type}_k${k}_jobs_${job_name}/job$job/log      # Make a directory to hold the log files (No rm command because the last one got it)

      dir=$(pwd)/${job_type}_k${k}_jobs_${job_name}/job$job     # The directory that the job executes in
      echo "Initialdir      = $dir" >>$file       # Set the different folders for job executions
      echo "Arguments       = \"${job} ${k} ${itype} \"" >>$file  # Set the different arguments (random seeds)
      echo "queue" >>$file                        # The job is good to go and condor will queue it

      # Now for the fun job of using a bash script to write a bash script
      echo "#!/bin/bash">>${job_type}_k${k}_jobs_${job_name}/job$job/run.sh          # Set up the bash file which will run your code
      echo "source $codedir/PythiaEnvSetup.sh">>${job_type}_k${k}_jobs_${job_name}/job$job/run.sh          # Source your setup file to make sure all requisite libraries are available
      echo "./$code \${1} \${2} \${3}">>${job_type}_k${k}_jobs_${job_name}/job$job/run.sh        # Execute the code with the argument condor provides

      echo "mv $output $outputpath/k${k}/${job_type}/pytree_\${1}.root">>${job_type}_k${k}_jobs_${job_name}/job$job/run.sh # Move the output tree to the right spot
      # echo "Done making k${k} ${job_type} job $job"
      chmod a+x ${job_type}_k${k}_jobs_${job_name}/job$job/run.sh     # Make the script we just wrote executable
    done
  done
done

condor_submit $file  # All done! Let condor work its magic.

cd $subm_script_dir
