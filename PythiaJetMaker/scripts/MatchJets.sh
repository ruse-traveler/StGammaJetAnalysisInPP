#!/bin/bash
# 'MatchJets.sh'
#
# Use this to run 'MatchJets.C' in batch mode

# filepaths
cPath="\"pythiaCheckP.r03a02rm1chrg.d27m1y2017.root\""
fPath="\"pythiaCheckP.r03a02rm1full.d27m1y2017.root\""
oPath="\"pythiaMatchP.r03a02rm1.d28m1y2017.root\""

root -b -q MatchJets.C"($cPath, $fPath, $oPath, true)"
