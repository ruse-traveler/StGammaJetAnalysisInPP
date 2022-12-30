#!/bin/bash
# 'CreateWeirdList.sh'
# Derek Anderson
#
# Use this to run 'CreateWeirdList.C'.
# Pass a list of files with "weird"
# events in them, and this will stream
# the file name, run ID, and event ID
# to a text file.


# input / output
list=$1
out="test.d9m4y2017.txt"


# make sure argument was passed
if [ -z $list ]; then
  printf "HEY! Please pass an argument!\n"
  exit
fi


# if so, read in list
for file in $(cat $list); do
  root -b -q "CreateWeirdList.C(\"$file\", \"$out\", false)"
done
