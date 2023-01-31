#!/bin/bash
# 'AutomateSubmission.sh'
# Derek Anderson
# 01.30.2023
#
# Use this to generate the necessary
# files and directories for running
# StJetTreeMcMaker on the grid
#
# NOTE: currently set-up for run14
# AuAu embedding sample.
#
# NOTE: $dir_to_process sets the
# directories to generate jobs
# for.  Should have a format
# similar to:
#   <path>/prefix_*_suffix/<rest_of_path>
# (The '*' is necessary.)
#
# NOTE: $sleep_duration sets how long
# b/n submitting waves of jobs

dir_to_process="/star/data105/embedding/AuAu_200_production_2014/Pythia6_pt9_11_*_20192901/P18ih.SL18h/2014"
sleep_duration="2s"

(( nDir=0 ))
for dir in $dir_to_process; do
  ./SubmitMatching.sh $dir $nDir
  sleep $sleep_duration
  (( nDir++ ))
done

# end -------------------------------------------------------------------------
