#!/bin/bash
# 'CombineWeirdLists.sh'
# Derek Anderson
#
# Have multiple lists of files with
# weird events?  Use this to combine
# them...


# lists to combine
declare -a lists
lists[0]="weirder.1digit.d6m4y2017.list"
lists[1]="weirder.2digit.d6m4y2017.list"
lists[2]="weirder.3digit.d6m4y2017.list"
lists[3]="weirder.4digit.d6m4y2017.list"

# parameters
out="weirdest.d6m4y2017.list"
pref="/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/tree.C27C92D9FEE17EE88FCFF22D7ECE052B_"
suff=".root"
nPre=${#pref}
nSuf=${#suff}


# organize file names
printf "\n  Beginning combining script...\n"
printf "    Organizing file names...\n"

nLine=0
iDig1=0
iDig2=0
iDig3=0
iDig4=0
declare -a digs1
declare -a digs2
declare -a digs3
declare -a digs4
for list in ${lists[@]}; do
  for line in $(cat $list); do
    nTot=${#line}
    nNot=$(expr $nPre + $nSuf)
    nDig=$(expr $nTot - $nNot)
    if [ $nDig == "1" ]; then
      digs1[iDig1]=$(echo $line)
      iDig1=$(( $iDig1 + 1 ))
    fi
    if [ $nDig == "2" ]; then
      digs2[iDig2]=$(echo $line)
      iDig2=$(( $iDig2 + 1 ))
    fi
    if [ $nDig == "3" ]; then
      digs3[iDig3]=$(echo $line)
      iDig3=$(( $iDig3 + 1 ))
    fi
    if [ $nDig == "4" ]; then
      digs4[iDig4]=$(echo $line)
      iDig4=$(( $iDig4 + 1 ))
    fi
    nLine=$(( $nLine + 1 ))
  done
done


# combine lists
printf "    Done organizing: $nLine lines sorted...\n"
printf "    Combining lists...\n"

nFile=0
iDig1=0
iDig2=0
iDig3=0
iDig4=0
declare -a files
# 1 digit files
for file in ${digs1[@]}; do
  if [ $iDig1 == "0" ]; then
    files[nFile]=$(echo $file)
    nFile=$(( $nFile + 1 ))
  else
    # check if file has already been added
    isSame=0
    for element in ${files[@]}; do
      if [ $file == $element ]; then
        isSame=1
        break
      fi
    done
    # if not, add to list
    if [ $isSame == "0" ]; then
      files[nFile]=$(echo $file)
      nFile=$(( $nFile + 1 ))
    fi
  fi
  iDig1=$(( $iDig1 + 1 ))
done
# 2 digit files
for file in ${digs2[@]}; do
  if [ $iDig2 == "0" ]; then
    files[nFile]=$(echo $file)
    nFile=$(( $nFile + 1 ))
  else
    # check if file has already been added
    isSame=0
    for element in ${files[@]}; do
      if [ $file == $element ]; then
        isSame=1
        break
      fi
    done
    # if not, add to list
    if [ $isSame == "0" ]; then
      files[nFile]=$(echo $file)
      nFile=$(( $nFile + 1 ))
    fi
  fi
  iDig2=$(( $iDig2 + 1 ))
done
# 3 digit files
for file in ${digs3[@]}; do
  if [ $iDig3 == "0" ]; then
    files[nFile]=$(echo $file)
    nFile=$(( $nFile + 1 ))
  else
    # check if file has already been added
    isSame=0
    for element in ${files[@]}; do
      if [ $file == $element ]; then
        isSame=1
        break
      fi
    done
    # if not, add to list
    if [ $isSame == "0" ]; then
      files[nFile]=$(echo $file)
      nFile=$(( $nFile + 1 ))
    fi
  fi
  iDig3=$(( $iDig3 + 1 ))
done
# 4 digit files
for file in ${digs4[@]}; do
  if [ $iDig4 == "0" ]; then
    files[nFile]=$(echo $file)
    nFile=$(( $nFile + 1 ))
  else
    # check if file has already been added
    isSame=0
    for element in ${files[@]}; do
      if [ $file == $element ]; then
        isSame=1
        break
      fi
    done
    # if not, add to list
    if [ $isSame == "0" ]; then
      files[nFile]=$(echo $file)
      nFile=$(( $nFile + 1 ))
    fi
  fi
  iDig4=$(( $iDig4 + 1 ))
done


# stream to output
printf "    Done combining lists: $nFile files added...\n"
printf "    Stream to output...\n"

for file in ${files[@]}; do
  printf "$file\n" >> $out
done


# delete arrays
unset lists
unset files
unset digs1
unset digs2
unset digs3
unset digs4

printf "  Combining script finished!\n\n"
