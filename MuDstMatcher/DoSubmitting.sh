#!/bin/bash
# 'DoSubmitting.sh'
#
# This script automates the submission process.  It
# will run the submission script for each partonic
# pT bin in staggered intervals.

# script names
declare -a scripts
scripts[0]="SubmitMatching9.sh"
scripts[1]="SubmitMatching11.sh"
scripts[2]="SubmitMatching15.sh"
scripts[3]="SubmitMatching25.sh"

# log names
declare -a logs
logs[0]="pt9.submit.log"
logs[1]="pt11.submit.log"
logs[2]="pt15.submit.log"
logs[3]="pt25.submit.log"


# begin script
printf "\nBeginning submission process!\n"

# delete old logs (if necesssary)
printf "  Checking for old logs.\n"
for log in ${logs[@]}; do
  if [ -f $log ]; then
    printf "WARNING: log '$log' exists! Deleting...\n"
    rm $log
  fi
done

# run scripts
(( nLog=0 ))
printf "  Running scripts.\n"
for script in ${scripts[@]}; do
  printf "    Starting script '$script'...\n"
  nohup ./$script >& ${logs[$nLog]} &
  printf "    Script '$script' finished! Going to sleep...\n"
  sleep 12h
  (( nLog++ ))
done
printf "  All scripts finished!\n"


# delete arrays
unset scripts
unset logs

# announce end
printf "Submission done!\n\n"
