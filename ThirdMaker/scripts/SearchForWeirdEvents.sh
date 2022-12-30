#!/bin/bash
# 'SearchForWeirdEvents.sh'
# Derek Anderson
#
# Use this to comb through the ThirdJetMaker
# output for weird events.


# i/o stuff
corpus="./tree.FB6ABC583E0B772ADF3EAAA8D675C833_*.out"
tarPre="   ----> variable "
tarSuf=" disagrees! <----"
outPre="tree.new.r10145070.var"
outSuf=".d17m4y2017.list"


for i in `seq 0 18`; do

  target=$tarPre$i$tarSuf
  output=$outPre$i$outSuf
  grep "$target" $corpus > $output

done
