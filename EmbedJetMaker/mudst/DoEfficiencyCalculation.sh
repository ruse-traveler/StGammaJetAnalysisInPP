#!/bin/bash
# 'DoEfficiencyCalculation.sh'
# Derek Anderson
# 02.02.2020
#
# Use this to run 'DoEfficiencyCalculation.C'
# in batch mode over a set of files


# declare arrays
declare -a input
declare -a output

# input files
input[0]="\"../../MuDstMatching/output/merged/pt5ff.matchWithMc.root\""
input[1]="\"../../MuDstMatching/output/merged/pt7ff.matchWithMc.root\""
input[2]="\"../../MuDstMatching/output/merged/pt9ff.matchWithMc.root\""
input[3]="\"../../MuDstMatching/output/merged/pt11ff.matchWithMc.root\""
input[4]="\"../../MuDstMatching/output/merged/pt15ff.matchWithMc.root\""
input[5]="\"../../MuDstMatching/output/merged/pt25ff.matchWithMc.root\""
input[6]="\"../../MuDstMatching/output/merged/pt35ff.matchWithMc.root\""
input[7]="\"../../MuDstMatching/output/merged/pt4rff.matchWithMc.root\""
input[8]="\"../../MuDstMatching/output/merged/pt5rff.matchWithMc.root\""
input[9]="\"../../MuDstMatching/output/merged/pt7rff.matchWithMc.root\""
input[10]="\"../../MuDstMatching/output/merged/pt9rff.matchWithMc.root\""
input[11]="\"../../MuDstMatching/output/merged/pt11rff.matchWithMc.root\""
input[12]="\"../../MuDstMatching/output/merged/pt15rff.matchWithMc.root\""
input[13]="\"../../MuDstMatching/output/merged/pt25rff.matchWithMc.root\""
input[14]="\"../../MuDstMatching/output/merged/pt35rff.matchWithMc.root\""

# output files
output[0]="\"pp200r9pt5ff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""
output[1]="\"pp200r9pt7ff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""
output[2]="\"pp200r9pt9ff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""
output[3]="\"pp200r9pt11ff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""
output[4]="\"pp200r9pt15ff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""
output[5]="\"pp200r9pt25ff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""
output[6]="\"pp200r9pt35ff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""
output[7]="\"pp200r9pt4rff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""
output[8]="\"pp200r9pt5rff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""
output[9]="\"pp200r9pt7rff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""
output[10]="\"pp200r9pt9rff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""
output[11]="\"pp200r9pt11rff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""
output[12]="\"pp200r9pt15rff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""
output[13]="\"pp200r9pt25rff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""
output[14]="\"pp200r9pt35rff.trackPseudoEff_withOnlyMatchVtxs_applyingOnlyPtGT0ForMatching.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root\""


# run macro
(( iFile=0 ))
for file in ${input[@]}; do
  root -b -q "DoEfficiencyCalculation.C($file, ${output[iFile]}, true)"
  (( iFile++ ))
done

# delete arrays
unset input
unset output
