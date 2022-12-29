#!/bin/bash
# 'GenerateDir.sh'
# Derek Anderson
#
# This generates directories for the output of 'SubmitMatching.sh'.

# arguments
output=$1
cwd=$PWD


# create general output dir.
tiptop=$cwd"/output"
if [ ! -d $tiptop ]; then
  printf "  GenerateDir.sh: Creating directory '$tiptop'\n"
  mkdir $tiptop
fi

# create specific output dir.
if [ ! -d $output ]; then
  printf "  GenerateDir.sh: Creating directory '$output'\n"
  mkdir $output
fi
