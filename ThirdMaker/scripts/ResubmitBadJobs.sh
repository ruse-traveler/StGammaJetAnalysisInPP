#!/bin/bash
# 'ResubmitBadJobs.sh'
# Derek Anderson
#
# Use this to fix output files with no keys. Run
# 'CheckTrees.C' first to generate a list of bad
# files.  Pass that list as an argument to this
# script.
#
# NOTE: that this needs the original '.dataset'
# file generated during the initial submission
# process in the same directory.


# input / output
list=$1
sess="/star/u/dmawxc/JetRecon_pp/FullJetTree/submit/pp200r9.d13m3y2017.session.xml"
pref="/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/pp200r9.d13m3y2017.job"
suff=".root"
nPre=${#pref}
nSuf=${#suff}
nNam=$(expr $nPre + $nSuf)


# make sure argument was passed
if [ -z $list ]; then
  printf "HEY! Please pass an argument!\n"
  exit
fi


# if so, read in list
printf "\nReading in '$list'...\n"

for file in $(cat $list); do

  # calculate job no.
  nTot=${#file}
  nNum=$(expr $nTot - $nNam)
  iJob=$(echo ${file:$nPre:$nNum})

  # submit job
  printf "  Submitting job '$iJob'...\n"
  star-submit -r $iJob $sess

done

printf "Submission script finished!\n\n"
