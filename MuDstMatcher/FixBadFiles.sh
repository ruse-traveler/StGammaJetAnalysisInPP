#!/bin/bash
# 'FixBadFiles.sh'
#
# Use this to fix output files with no keys.
# Run 'CheckEntries.C' first (in root), which
# will generate a list of bad files.  Pass that
# list as an argument to this script.


# input / output
ver="SL14a"
list=$1
lPath="./logs"
input="/projecta/projectdirs/starprod/embedding/production2009_200GeV/Jet_pp200_2009.elz17/SL11d_embed"
output="$PWD/output/pt15_25"
prefix="pt15_25_"
  
# various parameters
dStart=${#output}
pStart=${#prefix}
fStart=$(expr $dStart + 1)
rStart=$(expr $fStart + $pStart)
oldEnd="thirdmaker.root"
badEnd="badmaker.root"
dstEnd="MuDst.root"
newEnd="thirdmaker.root"
logEnd="fix.log"


printf "\nReading in '$list'...\n"

for file in $(cat $list); do

  # determine run and old filename
  run=$(echo ${file:$rStart:8})
  old=$(echo ${file:$fStart})
  bad=$(echo ${old/$oldEnd/$badEnd})
  dst=$(echo ${old/$oldEnd/$dstEnd})
  new=$(echo ${old/$oldEnd/$newEnd})
  log=$(echo ${old/$oldEnd/$logEnd})

  # create new filenames and arguments
  iPath=$input"/"$run
  oPath=$output
  oldPath=$oPath"/"$old
  badPath=$oPath"/"$bad
  dstPath=$iPath"/"$dst
  newPath=$oPath"/"$new
  logPath=$lPath"/"$log

  # generate new file
  if [ -a $oldPath ]; then
    mv $oldpath $badPath
    starver $ver
    root4star -b -q doFullJetTree.C\(2000000,\"$newPath\",\"$dstPath\"\) > $logPath;
  fi

  printf "  Regenerated file '$file'\n"

done  # file loop


printf "Finished!\n\n"
