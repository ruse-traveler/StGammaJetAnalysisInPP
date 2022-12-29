#!/bin/bash
# 'CalculateTrackEfficiency.sh'
# Derek Anderson
# 03.06.2019
#
# Use this to run 'CalculateTrackEfficiency.C'
# in batch mode


# declare arrays
declare -a input
declare -a output

# input files
input[0]="\"../../MuDstMatching/output/merged/pt4rff.matchWithMc.root\""
input[1]="\"../../MuDstMatching/output/merged/pt5rff.matchWithMc.root\""
input[2]="\"../../MuDstMatching/output/merged/pt7rff.matchWithMc.root\""
input[3]="\"../../MuDstMatching/output/merged/pt9rff.matchWithMc.root\""
input[4]="\"../../MuDstMatching/output/merged/pt11rff.matchWithMc.root\""
input[5]="\"../../MuDstMatching/output/merged/pt15rff.matchWithMc.root\""
input[6]="\"../../MuDstMatching/output/merged/pt25rff.matchWithMc.root\""
input[7]="\"../../MuDstMatching/output/merged/pt35rff.matchWithMc.root\""
input[8]="\"../../MuDstMatching/output/merged/pt5ff.matchWithMc.root\""
input[9]="\"../../MuDstMatching/output/merged/pt7ff.matchWithMc.root\""
input[10]="\"../../MuDstMatching/output/merged/pt9ff.matchWithMc.root\""
input[11]="\"../../MuDstMatching/output/merged/pt11ff.matchWithMc.root\""
input[12]="\"../../MuDstMatching/output/merged/pt15ff.matchWithMc.root\""
input[13]="\"../../MuDstMatching/output/merged/pt25ff.matchWithMc.root\""
input[14]="\"../../MuDstMatching/output/merged/pt35ff.matchWithMc.root\""

# output files
output[0]="\"pp200r9pt4rff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""
output[1]="\"pp200r9pt5rff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""
output[2]="\"pp200r9pt7rff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""
output[3]="\"pp200r9pt9rff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""
output[4]="\"pp200r9pt11rff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""
output[5]="\"pp200r9pt15rff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""
output[6]="\"pp200r9pt25rff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""
output[7]="\"pp200r9pt35rff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""
output[8]="\"pp200r9pt5ff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""
output[9]="\"pp200r9pt7ff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""
output[10]="\"pp200r9pt9ff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""
output[11]="\"pp200r9pt11ff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""
output[12]="\"pp200r9pt15ff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""
output[13]="\"pp200r9pt25ff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""
output[14]="\"pp200r9pt35ff.idEnd0ptMc0.et1220pt0230dca3vz55pi0.d10m5y2020.root\""

# run macro
(( iFile=0 ))
for file in ${input[@]}; do
  root -b -q "CalculateTrackEfficiency.C(${input[iFile]}, ${output[iFile]}, true)"
  (( iFile++ ))
done

# delete arrays
unset input
unset output
