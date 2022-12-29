#!/bin/bash
# 'SubmitJets.sh'
# Derek Anderson
#
# Use this to submit lots of jets!

cwd=$PWD
sub=$cwd"/submit"
sim="false"
ver="pro"
Mac=("MakeFemtoDst3.C" "MakeFemtoDst5.C")
Lvl=("p" "d")
nEvts=5302990
nFile=6
nEvtsPerFile=1000000

printf "Submitting $nEvts events...\n"

for mac in ${Mac[@]}; do
  for lvl in ${Lvl[@]}; do
    for i in `seq 0 5`; do

     # determine events
     j=$((  $i + 1 ))
     iStart=$(( $i*$nEvtsPerFile  ))
     iStop=$(( $j*$nEvtsPerFile - 1  ))
     if [ "$nEvts" -lt "$iStop" ]; then
       iStop=$nEvts
     fi

     # make output name
     if [ $mac = "MakeFemtoDst3.C" ]; then
       job="Pythia23"$lvl".r03a02rm1.gCharged"$i
       out=$cwd"/"$job".root"
     else
       job="Pythia23"$lvl".r05a065rm1.gCharged"$i
       out=$cwd"/"$job".root"
     fi

     # make argument
     if [ $lvl = "p" ]; then
       arg="'($iStart,$iStop,true,false,\"$out\")'"
     else
       arg="'($iStart,$iStop,false,true,\"$out\")'"
     fi

     cp $cwd"/"$mac $sub"/"
     ./GenerateXML.sh $sim $sub $ver $mac $arg $job
     mv $job".job.xml" $sub"/"
     cd $sub
     printf "  Submitting job '$job'...\n"
     printf "\n"
     star-submit $job".job.xml"
     printf "\n"
     cd $cwd

   done
 done
done

printf "Finished submitting!\n"
