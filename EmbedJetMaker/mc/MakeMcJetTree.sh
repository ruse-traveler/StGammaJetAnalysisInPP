#!/bin/bash
# 'MakeMuDstJetTree.sh'
# Derek Anderson
# 11.15.2019
#
# Use this to run 'MakeMuDstJetTree.C'
# in batch mode over a set of files


# declare arrays
declare -a pTpart
declare -a input
declare -a output

# partonic pT
pTpart[0]="5."
pTpart[1]="7."
pTpart[2]="9."
pTpart[3]="11."
pTpart[4]="15."
pTpart[5]="25."
pTpart[6]="35."
pTpart[7]="4."
pTpart[8]="5."
pTpart[9]="7."
pTpart[10]="9."
pTpart[11]="11."
pTpart[12]="15."
pTpart[13]="25."
pTpart[14]="35."

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
output[0]="\"pp200r9pt5ff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""
output[1]="\"pp200r9pt7ff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""
output[2]="\"pp200r9pt9ff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""
output[3]="\"pp200r9pt11ff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""
output[4]="\"pp200r9pt15ff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""
output[5]="\"pp200r9pt25ff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""
output[6]="\"pp200r9pt35ff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""
output[7]="\"pp200r9pt4rff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""
output[8]="\"pp200r9pt5rff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""
output[9]="\"pp200r9pt7rff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""
output[10]="\"pp200r9pt9rff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""
output[11]="\"pp200r9pt11rff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""
output[12]="\"pp200r9pt15rff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""
output[13]="\"pp200r9pt25rff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""
output[14]="\"pp200r9pt35rff.forUnfolding_Par.et920pt021Kvz55had.r05rm1chrg.d1m10y2020.root\""


# run macro
(( iFile=0 ))
for pT in ${pTpart[@]}; do
  root -b -q "MakeMcJetTree.C($pT, ${input[iFile]}, ${output[iFile]}, true)"
  (( iFile++ ))
done

# delete arrays
unset pTpart
unset input
unset output

# End -------------------------------------------------------------------------
