#!/bin/bash
# 'GenerateXML.sh'
# Derek Anderson
#
# This generates the job description file for star-submit.  Is called by the
# script 'SubmitMatching.sh'.

# arguments
xml=$1
sim=$2
num=$3
cwd=$4
ver=$5
mac=$6
arg=$7
lst=$8
log=$9


# warn if xml exists
name=$xml
if [ -f $name ]; then
  printf "GenerateXML.sh: WARNING! $name already exists.\n"
  printf "                Deleting and recreating...\n"
  rm $name
fi

# generate xml
touch $name
printf "<?xml version=\"1.0\" encoding=\"utf-8\" ?>\n" >> $name
printf "<job simulateSubmission=\"$sim\" fileListSyntax=\"paths\" maxFilesPerProcess=\"$num\">\n" >> $name
printf "\n" >> $name
printf "  <command>\n" >> $name
printf "    cd $cwd\n" >> $name
printf "    starver $ver\n" >> $name
printf "\n" >> $name
printf "    @ nFile=0\n" >> $name
printf '    while ( $nFile &lt; $INPUTFILECOUNT )'"\n" >> $name
printf "    eval set filename=\'\$INPUTFILE\'\$nFile\n" >> $name
printf "\n" >> $name
printf "    root4star -b -q $mac$arg\n" >> $name
printf "\n" >> $name
printf "    @ nFile++\n" >> $name
printf "    end\n" >> $name
printf "  </command>\n" >> $name
printf "\n" >> $name
printf "  <input URL=\"filelist:$cwd/$lst\" />\n" >> $name
printf "  <stdout URL=\"file:$log/\$JOBID.out\" />\n" >> $name
printf "  <stderr URL=\"file:$log/\$JOBID.err\" />\n" >> $name
printf "\n" >> $name
printf "</job>" >> $name
