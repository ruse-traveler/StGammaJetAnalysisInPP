#!/bin/bash
# 'ProcessWeirdList.sh'
# Derek Anderson
#
# Use this to process a list of 'weird' events
# (where some sort of event quantity differs
# between 'StMuEvent' and 'StEvent').  This
# will produce a list of files that contain
# 'weird' events.


# input / output
list=$1
path="/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/"
pref="tree.C27C92D9FEE17EE88FCFF22D7ECE052B_"
suff=".out"
ext=".root"
out="weirder.var0.d7m4y2017.list"
nPre=${#pref}
nSuf=${#suff}


# make sure argument was passed
if [ -z $list ]; then
  printf "HEY! Please pass an argument!\n"
  exit
fi


# if so, read in list
printf "\n  Reading in '$list'...\n"
printf "    Extracting files...\n"

iFile=0
declare -a files
for line in $(cat $list); do
  nName=${#line}
  files[$iFile]=$(echo ${line:0:$nName})
  iFile=$(( $iFile + 1 ))
done


# eliminate redundancies
printf "    Finished extracting files, reducing list...\n"

iFile=0
for file in ${files[@]}; do
  iPrev=$(( $iFile - 1 ))
  nName=${#files[iFile]}
  nBase=$(expr $nName - $nSuf)
  if [ "$iFile" == "0" ]; then
    base=$(echo ${files[iFile]:0:$nBase})
    name=$(echo $path$base$ext)
    printf "$name\n" >> $out
  else
    if [ "$file" != "${files[iPrev]}" ]; then
      base=$(echo ${files[iFile]:0:$nBase})
      name=$(echo $path$base$ext)
      printf "$name\n" >> $out
    fi
  fi
  iFile=$(( iFile + 1 ))
done
unset files

printf "  Processing script finished!\n\n"
