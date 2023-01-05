#!/bin/bash
# 'DoEfficiencyCalculation.sh'
# Derek Anderson
# 06.17.2020
#
# Use this to run 'DoTriggerEnergyScaleCalculation.C'
# in batch mode over a set of files

# declare arrays
declare -a input
declare -a output

# input files
#input[0]="\"./input/ThirdMaker/thirdMaker_oneEtaSite_merge.et630h09pi0.d25m5y2022.root\""
input[0]="\"./input/ThirdMaker/thirdMaker_merge.et650h09gam.d24m2y2021.root\""

# output files
output[0]="\"particleGun.forPwgComments_multiGamPerEvtNoVtx.et650gam.d18m12y2022.root\""

# tree names
parTree="\"McTracks\""
matTree="\"McTracksForMatching\""
detTree="\"Gfmtodst\""

# run macro
(( iFile=0 ))
for file in ${output[@]}; do
  root -b -q "DoTriggerEnergyScaleCalculation.C($file, ${input[iFile]}, ${input[iFile]}, $parTree, $matTree, $detTree, true)"
  (( iFile++ ))
done

# delete arrays
unset input
unset output

# End -------------------------------------------------------------------------
