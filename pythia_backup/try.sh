# job_name="default"
# if [ ${#1} -ne 0 ]
# then
#   job_name=${1}
# fi

if [ ${#1} -eq 0 ]
then
  echo "Must provide a non-empty string as argument!"
  echo "Suggestion: DATE_hardqcd_run# to avoid overwriting previous results"
  echo "e.g, 0314_hardqcd_run2"
  exit
fi

job_name=${1}
echo ${job_name}