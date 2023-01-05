#!/bin/bash
# 'MatchJets.sh'
#
# Use this to run 'MatchJets.C' in batch mode

# i/o filepaths
declare -a parFiles
declare -a detFiles
declare -a outFiles

# particle level input
parFiles[0]="\"../JetMaker/mc/output/March2020/pp200r9pt5ff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""
parFiles[1]="\"../JetMaker/mc/output/March2020/pp200r9pt7ff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""
parFiles[2]="\"../JetMaker/mc/output/March2020/pp200r9pt9ff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""
parFiles[3]="\"../JetMaker/mc/output/March2020/pp200r9pt11ff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""
parFiles[4]="\"../JetMaker/mc/output/March2020/pp200r9pt15ff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""
parFiles[5]="\"../JetMaker/mc/output/March2020/pp200r9pt25ff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""
parFiles[6]="\"../JetMaker/mc/output/March2020/pp200r9pt35ff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""
parFiles[7]="\"../JetMaker/mc/output/March2020/pp200r9pt4rff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""
parFiles[8]="\"../JetMaker/mc/output/March2020/pp200r9pt5rff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""
parFiles[9]="\"../JetMaker/mc/output/March2020/pp200r9pt7rff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""
parFiles[10]="\"../JetMaker/mc/output/March2020/pp200r9pt9rff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""
parFiles[11]="\"../JetMaker/mc/output/March2020/pp200r9pt11rff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""
parFiles[12]="\"../JetMaker/mc/output/March2020/pp200r9pt15rff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""
parFiles[13]="\"../JetMaker/mc/output/March2020/pp200r9pt25rff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""
parFiles[14]="\"../JetMaker/mc/output/March2020/pp200r9pt35rff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root\""

# detector level input
detFiles[0]="\"../JetMaker/mudst/output/March2020/pp200r9pt5ff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""
detFiles[1]="\"../JetMaker/mudst/output/March2020/pp200r9pt7ff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""
detFiles[2]="\"../JetMaker/mudst/output/March2020/pp200r9pt9ff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""
detFiles[3]="\"../JetMaker/mudst/output/March2020/pp200r9pt11ff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""
detFiles[4]="\"../JetMaker/mudst/output/March2020/pp200r9pt15ff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""
detFiles[5]="\"../JetMaker/mudst/output/March2020/pp200r9pt25ff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""
detFiles[6]="\"../JetMaker/mudst/output/March2020/pp200r9pt35ff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""
detFiles[7]="\"../JetMaker/mudst/output/March2020/pp200r9pt4rff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""
detFiles[8]="\"../JetMaker/mudst/output/March2020/pp200r9pt5rff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""
detFiles[9]="\"../JetMaker/mudst/output/March2020/pp200r9pt7rff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""
detFiles[10]="\"../JetMaker/mudst/output/March2020/pp200r9pt9rff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""
detFiles[11]="\"../JetMaker/mudst/output/March2020/pp200r9pt11rff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""
detFiles[12]="\"../JetMaker/mudst/output/March2020/pp200r9pt15rff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""
detFiles[13]="\"../JetMaker/mudst/output/March2020/pp200r9pt25rff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""
detFiles[14]="\"../JetMaker/mudst/output/March2020/pp200r9pt35rff.forUnfolding_Det.noTrigger.r02rm1chrg.d5m3y2020.root\""

# output
outFiles[0]="\"output/pp200r9pt5ff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""
outFiles[1]="\"output/pp200r9pt7ff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""
outFiles[2]="\"output/pp200r9pt9ff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""
outFiles[3]="\"output/pp200r9pt11ff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""
outFiles[4]="\"output/pp200r9pt15ff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""
outFiles[5]="\"output/pp200r9pt25ff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""
outFiles[6]="\"output/pp200r9pt35ff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""
outFiles[7]="\"output/pp200r9pt4rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""
outFiles[8]="\"output/pp200r9pt5rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""
outFiles[9]="\"output/pp200r9pt7rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""
outFiles[10]="\"output/pp200r9pt9rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""
outFiles[11]="\"output/pp200r9pt11rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""
outFiles[12]="\"output/pp200r9pt15rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""
outFiles[13]="\"output/pp200r9pt25rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""
outFiles[14]="\"output/pp200r9pt35rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root\""

# loop over files
(( iFile=0 ))
for pFile in ${parFiles[@]}; do

# run script
root -b -l <<EOF
  .L MatchJets.C++
  MatchJets($pFile, ${detFiles[$iFile]}, ${outFiles[$iFile]}, true)
  .q
EOF
(( iFile++ ))

# end file loop
done

# clean up files
rm *{.d,.so}

# delete arrays
unset parFiles
unset detFiles
unset outFiles
